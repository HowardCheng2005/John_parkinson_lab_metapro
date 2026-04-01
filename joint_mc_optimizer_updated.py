#!/usr/bin/env python3
import argparse
import csv
import math
import random
import re
from collections import Counter, defaultdict

READ_SUFFIX_RE = re.compile(r"(/[^|\s]+)$")


def normalize_read_key(read_id: str) -> str:
    return READ_SUFFIX_RE.sub("", str(read_id).strip())


def normalize_prediction_ec(prediction: str):
    if prediction is None:
        return None
    prediction = str(prediction).strip()
    if not prediction or prediction.lower() == "none":
        return None
    if prediction.startswith("EC:"):
        prediction = prediction[3:]
    return prediction


def softmax(scores, temp: float):
    if not scores:
        return []
    if temp <= 0:
        raise ValueError("Temperature must be > 0")
    max_score = max(scores)
    exps = [math.exp((score - max_score) / temp) for score in scores]
    total = sum(exps)
    if total <= 0:
        return [1.0 / len(scores)] * len(scores)
    return [value / total for value in exps]


def dedupe_best(candidates):
    best = {}
    for label, score in candidates:
        if label not in best or score > best[label]:
            best[label] = score
    return sorted(best.items(), key=lambda item: item[1], reverse=True)


def iter_deepec_txt_rows(path):
    """
    Parse a DeepEC plain-text result file.

    Expected data layout after the header:
        sequence_ID prediction score

    The parser is tolerant of:
    - spaces or tabs as delimiters
    - leading shell prompt / command lines before the header
    - blank lines
    """
    header_found = False

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            parts = line.split()
            if not header_found:
                if len(parts) >= 3 and parts[0] == "sequence_ID" and parts[1] == "prediction" and parts[2] == "score":
                    header_found = True
                continue

            if len(parts) < 3:
                continue

            yield parts[0], parts[1], parts[2]

    if not header_found:
        raise SystemExit(f"Could not find DeepEC header in {path}")


def load_metapathway(ec_path):
    pathway_ecs = set()

    with open(ec_path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue

            ec = parts[1].strip()
            if ec.lower() == "ec":
                continue
            pathway_ecs.add(ec)

    if not pathway_ecs:
        raise SystemExit(f"No EC values found in {ec_path}")

    return pathway_ecs


def load_deepec(deepec_path, pathway_ecs, topk_ec, ec_temp, keep_nonpath_ec=False):
    by_read = defaultdict(list)

    for sequence_id, prediction, score_text in iter_deepec_txt_rows(deepec_path):
        try:
            score = float(score_text)
        except ValueError:
            continue

        ec = normalize_prediction_ec(prediction)
        if ec is None:
            continue
        if ec not in pathway_ecs and not keep_nonpath_ec:
            continue

        read_key = normalize_read_key(sequence_id)
        by_read[read_key].append((ec, score))

    ec_candidates = {}
    for read_key, candidates in by_read.items():
        candidates = dedupe_best(candidates)
        if topk_ec > 0:
            candidates = candidates[:topk_ec]
        if not candidates:
            continue

        ec_labels = [ec for ec, _ in candidates]
        ec_scores = [score for _, score in candidates]
        ec_probs = softmax(ec_scores, ec_temp)
        ec_candidates[read_key] = (ec_labels, ec_probs)

    return ec_candidates


def load_centrifuger(centrifuger_path, topk_tax, tax_temp, tax_score_mode="margin"):
    by_read = defaultdict(list)

    with open(centrifuger_path, newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"readID", "taxID", "score"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise SystemExit(
                f"{centrifuger_path} must be a tab-delimited file with at least readID, taxID, and score columns"
            )

        for row in reader:
            read_id = row.get("readID")
            tax_id = row.get("taxID")
            score_text = row.get("score")
            second_text = row.get("2ndBestScore")

            if not read_id or not tax_id or score_text is None:
                continue

            try:
                score = float(score_text)
            except ValueError:
                continue

            score_used = score
            if tax_score_mode == "margin" and second_text not in (None, ""):
                try:
                    score_used = score - float(second_text)
                except ValueError:
                    score_used = score

            read_key = normalize_read_key(read_id)
            by_read[read_key].append((str(tax_id), score_used))

    tax_candidates = {}
    for read_key, candidates in by_read.items():
        candidates = dedupe_best(candidates)
        if topk_tax > 0:
            candidates = candidates[:topk_tax]
        if not candidates:
            continue

        tax_labels = [tax for tax, _ in candidates]
        tax_scores = [score for _, score in candidates]
        tax_probs = softmax(tax_scores, tax_temp)
        tax_candidates[read_key] = (tax_labels, tax_probs)

    return tax_candidates


def build_static_tax_abundance(tax_candidates):
    tax_mass = defaultdict(float)
    for tax_labels, tax_probs in tax_candidates.values():
        for tax, prob in zip(tax_labels, tax_probs):
            tax_mass[tax] += prob

    total_mass = sum(tax_mass.values())
    if total_mass <= 0:
        return {}

    return {tax: mass / total_mass for tax, mass in tax_mass.items()}


def metapathway_gain(tax_id, ec, seen_ecs, metapathway_ecs):
    if ec not in metapathway_ecs:
        return 0.0, 0.0

    covered = seen_ecs[tax_id]
    if ec in covered:
        return 0.0, 0.0

    total_ecs = (len(metapathway_ecs))
    gain = 1.0 / total_ecs
    repr_bonus = (len(covered) / total_ecs) * gain
    return gain, repr_bonus


def parse_args():
    parser = argparse.ArgumentParser(
        description="Monte Carlo joint taxon-EC optimizer using one KEGG metapathway"
    )
    parser.add_argument("-c", "--centrifuger", required=True,
                        help="Centrifuger output TSV with readID, taxID, score, and optionally 2ndBestScore")
    parser.add_argument("-d", "--deepec", required=True,
                        help="DeepEC plain-text result file with header: sequence_ID prediction score")
    parser.add_argument("-e", "--ec_tsv", required=True,
                        help="Two-column KEGG file: module<TAB>ec. All ECs are pooled into one metapathway")
    parser.add_argument("-o", "--outprefix", required=True,
                        help="Output prefix")

    parser.add_argument("--topk_tax", type=int, default=5,
                        help="Keep top-k taxon candidates per read. Use 0 or less to keep all")
    parser.add_argument("--topk_ec", type=int, default=5,
                        help="Keep top-k EC candidates per read. Use 0 or less to keep all")
    parser.add_argument("--tax_temp", type=float, default=20000.0,
                        help="Softmax temperature for taxon scores")
    parser.add_argument("--ec_temp", type=float, default=0.15,
                        help="Softmax temperature for EC scores")
    parser.add_argument("--tax_score_mode", choices=["raw", "margin"], default="margin",
                        help="Use raw tax score or score minus 2ndBestScore")

    parser.add_argument("-n", "--iters", type=int, default=200,
                        help="Number of Monte Carlo iterations")
    parser.add_argument("--seed", type=int, default=1,
                        help="Random seed")

    parser.add_argument("--lam_path", type=float, default=10.0,
                        help="Weight for adding a new EC to a taxon metapathway")
    parser.add_argument("--lam_repr", type=float, default=8.0,
                        help="Extra weight for extending an already represented metapathway")
    parser.add_argument("--lam_abund", type=float, default=2.0,
                        help="Weight for the static taxon abundance prior")

    parser.add_argument("--keep_nonpath_ec", action="store_true",
                        help="Keep DeepEC ECs not present in the KEGG EC file. They get no pathway bonus")

    args = parser.parse_args()

    if args.iters <= 0:
        parser.error("--iters must be > 0")
    if args.tax_temp <= 0:
        parser.error("--tax_temp must be > 0")
    if args.ec_temp <= 0:
        parser.error("--ec_temp must be > 0")

    return args


def main():
    args = parse_args()
    random.seed(args.seed)

    metapathway_ecs = load_metapathway(args.ec_tsv)
    ec_candidates = load_deepec(
        args.deepec,
        pathway_ecs=metapathway_ecs,
        topk_ec=args.topk_ec,
        ec_temp=args.ec_temp,
        keep_nonpath_ec=args.keep_nonpath_ec,
    )
    tax_candidates = load_centrifuger(
        args.centrifuger,
        topk_tax=args.topk_tax,
        tax_temp=args.tax_temp,
        tax_score_mode=args.tax_score_mode,
    )

    shared_reads = sorted(set(ec_candidates) & set(tax_candidates))
    if not shared_reads:
        raise SystemExit("No overlapping reads between DeepEC and Centrifuger after read-ID normalization")

    reads = []
    for read_key in shared_reads:
        ec_labels, ec_probs = ec_candidates[read_key]
        tax_labels, tax_probs = tax_candidates[read_key]
        if ec_labels and tax_labels:
            reads.append((read_key, tax_labels, tax_probs, ec_labels, ec_probs))

    if not reads:
        raise SystemExit("No reads left after filtering. Each retained read must have at least one taxon and one EC candidate")

    tax_abundance = build_static_tax_abundance(tax_candidates)

    joint_counts = defaultdict(Counter)
    tax_counts = defaultdict(Counter)
    ec_counts = defaultdict(Counter)
    tax_metapath_sum = defaultdict(float)
    tax_seen_iters = Counter()

    for _ in range(args.iters):
        shuffled_reads = list(reads)
        random.shuffle(shuffled_reads)
        seen_ecs = defaultdict(set)

        for read_key, tax_labels, tax_probs, ec_labels, ec_probs in shuffled_reads:
            pair_labels = []
            pair_weights = []

            for tax, tax_prior in zip(tax_labels, tax_probs):
                abundance_bonus = tax_abundance.get(tax, 0.0)
                for ec, ec_prior in zip(ec_labels, ec_probs):
                    path_gain, repr_bonus = metapathway_gain(tax, ec, seen_ecs, metapathway_ecs)
                    log_weight = math.log(max(tax_prior, 1e-300)) + math.log(max(ec_prior, 1e-300))
                    log_weight += args.lam_path * path_gain
                    log_weight += args.lam_repr * repr_bonus
                    log_weight += args.lam_abund * abundance_bonus

                    pair_labels.append((tax, ec))
                    pair_weights.append(math.exp(log_weight))

            chosen_tax, chosen_ec = random.choices(pair_labels, weights=pair_weights, k=1)[0]
            joint_counts[read_key][(chosen_tax, chosen_ec)] += 1
            tax_counts[read_key][chosen_tax] += 1
            ec_counts[read_key][chosen_ec] += 1

            if chosen_ec in metapathway_ecs:
                seen_ecs[chosen_tax].add(chosen_ec)

        for tax, covered_ecs in seen_ecs.items():
            tax_seen_iters[tax] += 1
            completeness = len(covered_ecs) / float(len(metapathway_ecs))
            tax_metapath_sum[tax] += completeness

    enzyme_taxon_counts = defaultdict(Counter)
    all_taxa = set()
    all_ecs = set()

    with open(f"{args.outprefix}.read_joint_map.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["readID", "MAP_taxID", "MAP_EC", "joint_prob", "tax_marginal_prob", "ec_marginal_prob"])
        for read_key, *_ in reads:
            (best_tax, best_ec), best_joint_count = joint_counts[read_key].most_common(1)[0]

            enzyme_taxon_counts[best_ec][best_tax] += 1
            all_taxa.add(best_tax)
            all_ecs.add(best_ec)

            writer.writerow([
                read_key,
                best_tax,
                best_ec,
                f"{best_joint_count / args.iters:.6f}",
                f"{tax_counts[read_key][best_tax] / args.iters:.6f}",
                f"{ec_counts[read_key][best_ec] / args.iters:.6f}",
            ])

    with open(f"{args.outprefix}.enzyme_taxon_matrix.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        sorted_taxa = sorted(all_taxa)
        writer.writerow(["EC"] + sorted_taxa)
        for ec in sorted(all_ecs):
            writer.writerow([ec] + [enzyme_taxon_counts[ec].get(tax, 0) for tax in sorted_taxa])

    with open(f"{args.outprefix}.read_joint_all.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["readID", "taxID", "EC", "posterior_prob"])
        for read_key in sorted(joint_counts):
            for (tax, ec), count in joint_counts[read_key].most_common():
                writer.writerow([read_key, tax, ec, f"{count / args.iters:.6f}"])

    with open(f"{args.outprefix}.taxon_metapathway_completeness.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["taxID", "mean_metapathway_completeness", "metapathway_size", "iterations_observed"])
        for tax, total in sorted(tax_metapath_sum.items(), key=lambda item: item[1], reverse=True):
            writer.writerow([
                tax,
                f"{total / args.iters:.6f}",
                len(metapathway_ecs),
                tax_seen_iters[tax],
            ])

    with open(f"{args.outprefix}.candidate_priors.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["readID", "candidate_type", "label", "prior_prob"])
        for read_key in sorted(tax_candidates):
            tax_labels, tax_probs = tax_candidates[read_key]
            for tax, prob in zip(tax_labels, tax_probs):
                writer.writerow([read_key, "tax", tax, f"{prob:.6f}"])
        for read_key in sorted(ec_candidates):
            ec_labels, ec_probs = ec_candidates[read_key]
            for ec, prob in zip(ec_labels, ec_probs):
                writer.writerow([read_key, "ec", ec, f"{prob:.6f}"])


if __name__ == "__main__":
    main()