#!/usr/bin/env python3
import argparse
import csv
import re
from collections import defaultdict

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


def iter_deepec_txt_rows(path):
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
    if not ec_path:
        return None

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


def load_best_deepec(deepec_path, pathway_ecs=None, keep_nonpath_ec=False):
    best_by_read = {}

    for sequence_id, prediction, score_text in iter_deepec_txt_rows(deepec_path):
        try:
            score = float(score_text)
        except ValueError:
            continue

        ec = normalize_prediction_ec(prediction)
        if ec is None:
            continue
        if pathway_ecs is not None and ec not in pathway_ecs and not keep_nonpath_ec:
            continue

        read_key = normalize_read_key(sequence_id)
        if read_key not in best_by_read or score > best_by_read[read_key][1]:
            best_by_read[read_key] = (ec, score)

    return best_by_read


def load_best_centrifuger(centrifuger_path, tax_score_mode="margin"):
    best_by_read = {}

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
            if read_key not in best_by_read or score_used > best_by_read[read_key][1]:
                best_by_read[read_key] = (str(tax_id), score_used)

    return best_by_read


def parse_args():
    parser = argparse.ArgumentParser(
        description="Deterministic top-hit taxon-EC baseline"
    )
    parser.add_argument("-c", "--centrifuger", required=True,
                        help="Centrifuger output TSV with readID, taxID, score, and optionally 2ndBestScore")
    parser.add_argument("-d", "--deepec", required=True,
                        help="DeepEC plain-text result file with header: sequence_ID prediction score")
    parser.add_argument("-o", "--outprefix", required=True,
                        help="Output prefix")
    parser.add_argument("-e", "--ec_tsv", default=None,
                        help="Optional two-column KEGG file: module<TAB>ec. If provided, completeness uses this metapathway")
    parser.add_argument("--tax_score_mode", choices=["raw", "margin"], default="margin",
                        help="Use raw tax score or score minus 2ndBestScore")
    parser.add_argument("--keep_nonpath_ec", action="store_true",
                        help="If --ec_tsv is provided, keep DeepEC ECs not present in the KEGG EC file")
    return parser.parse_args()


def main():
    args = parse_args()

    metapathway_ecs = load_metapathway(args.ec_tsv)
    best_ec_by_read = load_best_deepec(
        args.deepec,
        pathway_ecs=metapathway_ecs,
        keep_nonpath_ec=args.keep_nonpath_ec,
    )
    best_tax_by_read = load_best_centrifuger(
        args.centrifuger,
        tax_score_mode=args.tax_score_mode,
    )

    shared_reads = sorted(set(best_ec_by_read) & set(best_tax_by_read))
    if not shared_reads:
        raise SystemExit("No overlapping reads between DeepEC and Centrifuger after read-ID normalization")

    enzyme_taxon_counts = defaultdict(lambda: defaultdict(int))
    tax_to_ecs = defaultdict(set)
    all_taxa = set()
    all_ecs = set()

    for read_key in shared_reads:
        best_tax, _ = best_tax_by_read[read_key]
        best_ec, _ = best_ec_by_read[read_key]

        enzyme_taxon_counts[best_ec][best_tax] += 1
        tax_to_ecs[best_tax].add(best_ec)
        all_taxa.add(best_tax)
        all_ecs.add(best_ec)

    if metapathway_ecs is not None:
        completeness_universe = set(metapathway_ecs)
    else:
        completeness_universe = set(all_ecs)

    metapathway_size = len(completeness_universe)
    if metapathway_size == 0:
        raise SystemExit("No EC values available to summarize")

    with open(f"{args.outprefix}.read_joint_map.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["readID", "MAP_taxID", "MAP_EC", "joint_prob", "tax_marginal_prob", "ec_marginal_prob"])
        for read_key in shared_reads:
            best_tax, _ = best_tax_by_read[read_key]
            best_ec, _ = best_ec_by_read[read_key]
            writer.writerow([read_key, best_tax, best_ec, "1.000000", "1.000000", "1.000000"])

    with open(f"{args.outprefix}.enzyme_taxon_matrix.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        sorted_taxa = sorted(all_taxa)
        writer.writerow(["EC"] + sorted_taxa)
        for ec in sorted(all_ecs):
            writer.writerow([ec] + [enzyme_taxon_counts[ec].get(tax, 0) for tax in sorted_taxa])

    with open(f"{args.outprefix}.read_joint_all.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["readID", "taxID", "EC", "posterior_prob"])
        for read_key in shared_reads:
            best_tax, _ = best_tax_by_read[read_key]
            best_ec, _ = best_ec_by_read[read_key]
            writer.writerow([read_key, best_tax, best_ec, "1.000000"])

    with open(f"{args.outprefix}.taxon_metapathway_completeness.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["taxID", "mean_metapathway_completeness", "metapathway_size", "iterations_observed"])
        for tax, ecs in sorted(tax_to_ecs.items(), key=lambda item: len(item[1]), reverse=True):
            overlap = len(ecs & completeness_universe)
            writer.writerow([tax, f"{overlap / float(metapathway_size):.6f}", metapathway_size, 1])

    with open(f"{args.outprefix}.candidate_priors.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["readID", "candidate_type", "label", "prior_prob"])
        for read_key in shared_reads:
            best_tax, _ = best_tax_by_read[read_key]
            best_ec, _ = best_ec_by_read[read_key]
            writer.writerow([read_key, "tax", best_tax, "1.000000"])
            writer.writerow([read_key, "ec", best_ec, "1.000000"])


if __name__ == "__main__":
    main()