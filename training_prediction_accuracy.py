#!/usr/bin/env python3
import argparse
import csv
import re
from collections import Counter, defaultdict

EC_TAG_RE = re.compile(r"\|EC=([^|\s]+)")
READ_SUFFIX_RE = re.compile(r"(/[^|\s]+)$")


def normalize_read_id(read_id: str) -> str:
    return READ_SUFFIX_RE.sub("", str(read_id).strip())


def extract_true_ec_from_read_id(read_id: str):
    """
    Extract true EC from read ID like:
    000000001|WP_xxx|EC=1.2.1.11_0_0
    -> 1.2.1.11

    000000001|WP_xxx|EC=-1_0_0
    -> -1
    """
    m = EC_TAG_RE.search(str(read_id))
    if not m:
        return None
    return m.group(1).split("_", 1)[0]


def normalize_pred_ec(pred):
    if pred is None:
        return None
    pred = str(pred).strip()
    if not pred or pred.lower() == "none":
        return None
    if pred.startswith("EC:"):
        pred = pred[3:]
    return pred


def iter_deepec_rows(path):
    """
    Parse DeepEC plain-text output:
    sequence_ID prediction score
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


def detect_tsv_columns(fieldnames, user_read_col=None, user_pred_col=None):
    if fieldnames is None:
        raise SystemExit("Input TSV has no header")

    fields = set(fieldnames)

    read_candidates = ["readID", "sequence_ID", "read_id", "seqID", "seq_id"]
    pred_candidates = ["MAP_EC", "prediction", "predicted_EC", "EC_pred", "ec_pred", "EC"]

    read_col = user_read_col if user_read_col else next((c for c in read_candidates if c in fields), None)
    pred_col = user_pred_col if user_pred_col else next((c for c in pred_candidates if c in fields), None)

    if read_col is None:
        raise SystemExit(
            f"Could not find read ID column. Available columns: {', '.join(fieldnames)}"
        )
    if pred_col is None:
        raise SystemExit(
            f"Could not find predicted EC column. Available columns: {', '.join(fieldnames)}"
        )

    return read_col, pred_col


def iter_tsv_rows(path, read_col=None, pred_col=None):
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        read_col, pred_col = detect_tsv_columns(reader.fieldnames, read_col, pred_col)

        for row in reader:
            rid = row.get(read_col)
            pred = row.get(pred_col)
            if rid is None or pred is None:
                continue
            yield rid, pred


def load_predictions(path, fmt="auto", read_col=None, pred_col=None):
    """
    Returns list of tuples:
    (normalized_read_id, true_ec, predicted_ec)
    """
    rows = []

    if fmt == "deepec":
        iterator = iter_deepec_rows(path)
        for read_id, prediction, _score in iterator:
            read_id_norm = normalize_read_id(read_id)
            true_ec = extract_true_ec_from_read_id(read_id)
            pred_ec = normalize_pred_ec(prediction)
            rows.append((read_id_norm, true_ec, pred_ec))

    elif fmt == "tsv":
        iterator = iter_tsv_rows(path, read_col=read_col, pred_col=pred_col)
        for read_id, prediction in iterator:
            read_id_norm = normalize_read_id(read_id)
            true_ec = extract_true_ec_from_read_id(read_id)
            pred_ec = normalize_pred_ec(prediction)
            rows.append((read_id_norm, true_ec, pred_ec))

    else:
        # auto-detect
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            first_nonempty = None
            for raw in handle:
                line = raw.strip()
                if line:
                    first_nonempty = line
                    break

        if first_nonempty is None:
            raise SystemExit(f"Input file is empty: {path}")

        if first_nonempty.split()[:3] == ["sequence_ID", "prediction", "score"]:
            return load_predictions(path, fmt="deepec", read_col=read_col, pred_col=pred_col)
        else:
            return load_predictions(path, fmt="tsv", read_col=read_col, pred_col=pred_col)

    return rows


def main():
    ap = argparse.ArgumentParser(
        description="Compare predicted ECs against ground-truth EC tags embedded in read IDs"
    )
    ap.add_argument("-i", "--input", required=True,
                    help="Prediction file: DeepEC text output or TSV such as read_joint_map.tsv")
    ap.add_argument("-o", "--outprefix", required=True,
                    help="Output prefix")
    ap.add_argument("--format", choices=["auto", "deepec", "tsv"], default="auto",
                    help="Input format")
    ap.add_argument("--read_col", default=None,
                    help="Read ID column name for TSV input")
    ap.add_argument("--pred_col", default=None,
                    help="Predicted EC column name for TSV input")
    ap.add_argument("--skip_minus1", action="store_true",
                    help="Skip reads whose true EC tag is -1")
    ap.add_argument("--ignore_missing_pred", action="store_true",
                    help="Skip rows where predicted EC is missing/None instead of counting them as wrong")
    args = ap.parse_args()

    rows = load_predictions(
        args.input,
        fmt=args.format,
        read_col=args.read_col,
        pred_col=args.pred_col,
    )

    total = 0
    correct = 0
    incorrect = 0
    skipped_no_truth = 0
    skipped_minus1 = 0
    skipped_missing_pred = 0

    per_true_total = Counter()
    per_true_correct = Counter()
    mismatch_counts = Counter()

    mismatch_rows = []

    for read_id, true_ec, pred_ec in rows:
        if true_ec is None:
            skipped_no_truth += 1
            continue

        if args.skip_minus1 and true_ec == "-1":
            skipped_minus1 += 1
            continue

        if pred_ec is None and args.ignore_missing_pred:
            skipped_missing_pred += 1
            continue

        total += 1
        per_true_total[true_ec] += 1

        if pred_ec == true_ec:
            correct += 1
            per_true_correct[true_ec] += 1
        else:
            incorrect += 1
            mismatch_counts[(true_ec, pred_ec)] += 1
            mismatch_rows.append((read_id, true_ec, pred_ec if pred_ec is not None else "None"))

    accuracy = (correct / total) if total > 0 else 0.0

    with open(f"{args.outprefix}.summary.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["total_evaluated", total])
        writer.writerow(["correct", correct])
        writer.writerow(["incorrect", incorrect])
        writer.writerow(["accuracy", f"{accuracy:.6f}"])
        writer.writerow(["skipped_no_truth_tag", skipped_no_truth])
        writer.writerow(["skipped_minus1", skipped_minus1])
        writer.writerow(["skipped_missing_pred", skipped_missing_pred])

    with open(f"{args.outprefix}.per_ec.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["true_ec", "n", "correct", "accuracy"])
        for ec in sorted(per_true_total):
            n = per_true_total[ec]
            c = per_true_correct[ec]
            acc = c / n if n > 0 else 0.0
            writer.writerow([ec, n, c, f"{acc:.6f}"])

    with open(f"{args.outprefix}.mismatches.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["readID", "true_ec", "pred_ec"])
        for row in mismatch_rows:
            writer.writerow(row)

    with open(f"{args.outprefix}.mismatch_counts.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["true_ec", "pred_ec", "count"])
        for (true_ec, pred_ec), count in mismatch_counts.most_common():
            writer.writerow([true_ec, pred_ec if pred_ec is not None else "None", count])

    print(f"Total evaluated: {total}")
    print(f"Correct: {correct}")
    print(f"Incorrect: {incorrect}")
    print(f"Accuracy: {accuracy:.6f}")
    print(f"Skipped no truth tag: {skipped_no_truth}")
    print(f"Skipped true EC = -1: {skipped_minus1}")
    print(f"Skipped missing prediction: {skipped_missing_pred}")


if __name__ == "__main__":
    main()