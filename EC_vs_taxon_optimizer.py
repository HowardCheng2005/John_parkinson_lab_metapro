#!/usr/bin/env python3
import argparse, csv, random, math, re
from collections import defaultdict, Counter

EC_RE = re.compile(r"\|EC=([^|\s]+)")

def extract_ec(read_id: str):
    m = EC_RE.search(str(read_id))
    return None if not m else m.group(1).split("_", 1)[0]

def softmax(scores, temp):
    mx = max(scores)
    exps = [math.exp((s - mx) / temp) for s in scores]
    s = sum(exps)
    return [e / s for e in exps]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--centrifuger", required=True, help="centrifuger output TSV (with -k hits)")
    ap.add_argument("-e", "--ec_tsv", required=True, help="module-EC TSV: module<TAB>EC")
    ap.add_argument("-m", "--module", required=True, help="module id (e.g. M00001)")
    ap.add_argument("-k", type=int, default=5, help="top-k candidates per read (default 5)")
    ap.add_argument("-n", "--iters", type=int, default=200, help="MC iterations (default 200)")
    ap.add_argument("--temp", type=float, default=20000.0, help="softmax temperature for scores (default 20000)")
    ap.add_argument("--lam", type=float, default=10.0, help="pathway gain strength (default 10)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--skip_minus1", action="store_true", help="skip reads with EC=-1")
    ap.add_argument("-o", "--outprefix", required=True, help="output prefix")
    args = ap.parse_args()

    random.seed(args.seed)

    # 1) Load module EC set
    module_ecs = set()
    with open(args.ec_tsv, newline="") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            mod, ec = line.split("\t")[:2]
            if mod == args.module:
                module_ecs.add(ec)
    if not module_ecs:
        raise SystemExit(f"No ECs found for module {args.module} in {args.ec_tsv}")
    denom = float(len(module_ecs))  # completeness increment = 1/denom

    # 2) Load centrifuger candidates grouped by read
    by_read = defaultdict(list)  # readID -> list[(taxID, score)]
    read_ec = {}                 # readID -> EC

    with open(args.centrifuger, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rid = row["readID"]
            ec = extract_ec(rid)
            if ec is None:
                continue
            if args.skip_minus1 and ec == "-1":
                continue
            if ec not in module_ecs:
                continue

            tax = str(row["taxID"])
            score = int(float(row["score"]))
            read_ec[rid] = ec
            by_read[rid].append((tax, score))

    reads = []
    for rid, lst in by_read.items():
        lst.sort(key=lambda x: x[1], reverse=True)
        lst = lst[:args.k]
        taxes = [t for t, _ in lst]
        scores = [s for _, s in lst]
        priors = softmax(scores, args.temp)   # base probs from ranking
        reads.append((rid, read_ec[rid], taxes, priors))

    if not reads:
        raise SystemExit("No reads left after filtering to module ECs (check inputs).")

    # 3) Monte Carlo
    pick_counts = defaultdict(Counter)  # readID -> taxID -> count
    tax_score_sum = defaultdict(float)  # taxID -> sum completeness across iterations

    for _ in range(args.iters):
        random.shuffle(reads)

        seen = defaultdict(set)        # taxID -> set(EC already counted this iter)
        iter_score = defaultdict(float)

        for rid, ec, taxes, priors in reads:
            weights = []
            for t, p in zip(taxes, priors):
                gain = (1.0/denom) if ec not in seen[t] else 0.0
                weights.append(p * math.exp(args.lam * gain))

            chosen = random.choices(taxes, weights=weights, k=1)[0]
            pick_counts[rid][chosen] += 1

            if ec not in seen[chosen]:
                seen[chosen].add(ec)
                iter_score[chosen] += 1.0/denom

        for t, sc in iter_score.items():
            tax_score_sum[t] += sc

    # 4) Outputs
    with open(f"{args.outprefix}.read_map.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["readID", "EC", "MAP_taxID", "MAP_prob"])
        for rid, ec, _, _ in reads:
            best_tax, best_ct = pick_counts[rid].most_common(1)[0]
            w.writerow([rid, ec, best_tax, f"{best_ct/args.iters:.4f}"])

    with open(f"{args.outprefix}.taxon_mean_completeness.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["taxID", "mean_completeness"])
        for t, s in sorted(tax_score_sum.items(), key=lambda x: x[1], reverse=True):
            w.writerow([t, f"{s/args.iters:.6f}"])

if __name__ == "__main__":
    main()