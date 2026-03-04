#!/usr/bin/env python3
import argparse, csv, random, math, re
from collections import defaultdict, Counter

# Regex to extract EC from readID strings like: 000000001|WP_xxx|EC=1.2.1.11_0_0
EC_RE = re.compile(r"\|EC=([^|\s]+)")

def extract_ec(read_id: str):
    """
    Extract the EC number from a readID.
    Returns None if no EC tag is found.
    Example: "...|EC=1.2.1.11_0_0" -> "1.2.1.11"
             "...|EC=-1_0_0"       -> "-1"
    """
    m = EC_RE.search(str(read_id))
    return None if not m else m.group(1).split("_", 1)[0]

def softmax(scores, temp: float):
    """
    Convert a list of scores into probabilities using a stable softmax.
    Higher temp => more uniform probabilities.
    Lower temp  => more concentrated on the top score.
    """
    mx = max(scores)
    exps = [math.exp((s - mx) / temp) for s in scores]
    total = sum(exps)
    return [e / total for e in exps]

def main():
    ap = argparse.ArgumentParser()

    # Inputs
    ap.add_argument("-c", "--centrifuger", required=True,
                    help="Centrifuger output TSV (must include readID,taxID,score; use -k for multiple hits)")
    ap.add_argument("-e", "--EC", required=True,
                    help="Module-EC TSV with 2 columns: module<TAB>ec (ALL modules are used)")

    # MC + scoring settings
    ap.add_argument("-k", type=int, default=5,
                    help="Use top-k taxon candidates per read (default 5)")
    ap.add_argument("-n", "--iters", type=int, default=200,
                    help="Monte Carlo iterations (default 200)")
    ap.add_argument("--temp", type=float, default=20000.0,
                    help="Softmax temperature for score->prob conversion (default 20000)")
    ap.add_argument("--lam", type=float, default=10.0,
                    help="Strength of pathway-completeness bias (default 10)")
    ap.add_argument("--seed", type=int, default=1,
                    help="Random seed for reproducibility")

    # Filtering
    ap.add_argument("--skip_minus1", action="store_true",
                    help="Skip reads with EC=-1")

    # Output
    ap.add_argument("-o", "--outprefix", required=True,
                    help="Output prefix (no extension)")

    args = ap.parse_args()
    random.seed(args.seed)

    # ------------------------------------------------------------
    # 1) Load ALL ECs from the EC TSV into a set (no pandas)
    #    This defines the "pathway EC universe" we care about.
    # ------------------------------------------------------------
    all_ec = set()
    with open(args.EC, newline="") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Expect at least: module<TAB>ec
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            ec_val = parts[1].strip()
            if ec_val:
                all_ec.add(ec_val)

    if not all_ec:
        raise SystemExit(f"No ECs found in {args.EC}")

    # Each new EC covered contributes 1/|all_ec| to completeness
    denom = float(len(all_ec))

    # ------------------------------------------------------------
    # 2) Load centrifuger hits grouped by readID (pathway-only reads)
    #    - Only keep reads whose extracted EC is in all_ec
    #    - Store candidate (taxID, score) pairs per readID
    # ------------------------------------------------------------
    by_read = defaultdict(list)  # readID -> list of (taxID, score)
    read_ec = {}                 # readID -> extracted EC

    with open(args.centrifuger, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rid = row.get("readID")
            tax = row.get("taxID")
            score_str = row.get("score")

            if rid is None or tax is None or score_str is None:
                continue

            ec = extract_ec(rid)
            if ec is None:
                continue
            if args.skip_minus1 and ec == "-1":
                continue

            # PATHWAY-ONLY filter: keep only reads whose EC is in the EC set
            if ec not in all_ec:
                continue

            try:
                score = int(float(score_str))
            except Exception:
                continue

            read_ec[rid] = ec
            by_read[rid].append((str(tax), score))

    # ------------------------------------------------------------
    # 3) Convert each read's candidates into:
    #    (readID, EC, [taxIDs], [base_probs])
    #    where base_probs is softmax(score) over top-k hits.
    # ------------------------------------------------------------
    reads = []
    for rid, cand_list in by_read.items():
        # Sort candidates by descending score and keep top-k
        cand_list.sort(key=lambda x: x[1], reverse=True)
        cand_list = cand_list[:args.k]

        taxes = [t for t, _ in cand_list]
        scores = [s for _, s in cand_list]

        # Base probabilities from ranking/score only
        priors = softmax(scores, args.temp)

        reads.append((rid, read_ec[rid], taxes, priors))

    if not reads:
        raise SystemExit("No reads left after filtering to pathway ECs. Check EC list and readID EC tags.")

    # ------------------------------------------------------------
    # 4) Monte Carlo simulation
    #
    # For each iteration:
    #   - Shuffle reads (avoid order bias)
    #   - Track, for each taxon, which ECs have been "covered" this iteration (seen[taxon])
    #   - For each read, sample a taxon among its candidates with weight:
    #         weight = base_prob * exp(lam * gain)
    #     where gain = 1/denom if assigning this EC to this taxon adds a new EC (improves completeness),
    #           gain = 0 otherwise.
    #
    # We record:
    #   - pick_counts[readID][taxID] = how often each taxID was chosen for that read
    #   - tax_score_sum[taxID] = sum of completeness contributions across iterations
    # ------------------------------------------------------------
    pick_counts = defaultdict(Counter)     # readID -> Counter(taxID -> count)
    tax_score_sum = defaultdict(float)     # taxID -> sum completeness across iters

    for _ in range(args.iters):
        random.shuffle(reads)

        seen = defaultdict(set)           # taxID -> set of ECs counted this iteration
        iter_score = defaultdict(float)   # taxID -> completeness this iteration

        for rid, ec, taxes, priors in reads:
            weights = []
            for t, p in zip(taxes, priors):
                # Gain if ec is newly added to this taxon's covered EC set
                gain = (1.0 / denom) if ec not in seen[t] else 0.0
                # Bias the sampling toward choices that improve completeness
                weights.append(p * math.exp(args.lam * gain))

            chosen_tax = random.choices(taxes, weights=weights, k=1)[0]
            pick_counts[rid][chosen_tax] += 1

            # If new EC for that taxon, update completeness trackers
            if ec not in seen[chosen_tax]:
                seen[chosen_tax].add(ec)
                iter_score[chosen_tax] += 1.0 / denom

        # Add this iteration's completeness to totals
        for t, sc in iter_score.items():
            tax_score_sum[t] += sc

    # ------------------------------------------------------------
    # 5) Output results
    #
    # (A) Per read: MAP_taxID (most frequently chosen) + MAP_prob
    # (B) Per taxon: mean completeness across iterations
    # ------------------------------------------------------------
    with open(f"{args.outprefix}.read_map.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["readID", "EC", "MAP_taxID", "MAP_prob"])
        for rid, ec, _, _ in reads:
            best_tax, best_ct = pick_counts[rid].most_common(1)[0]
            w.writerow([rid, ec, best_tax, f"{best_ct/args.iters:.4f}"])

    with open(f"{args.outprefix}.taxon_mean_completeness.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["taxID", "mean_completeness"])
        for tax_id, total in sorted(tax_score_sum.items(), key=lambda x: x[1], reverse=True):
            w.writerow([tax_id, f"{total/args.iters:.6f}"])

if __name__ == "__main__":
    main()