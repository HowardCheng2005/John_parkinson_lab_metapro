import argparse
import pandas as pd
import re

EC_RE = re.compile(r"\|EC=([^|\s]+)")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--centrifuger", required=True, help="input centrifuger file")
    parser.add_argument("-e", "--EC", required=True, help="input EC file")
    parser.add_argument("-o", "--output", required=True, help="output file")
    args = parser.parse_args()

    taxons = pd.read_csv(args.centrifuger, sep="\t")
    ec = pd.read_csv(args.EC, sep="\t", header=None, names=["module", "ec"])

    all_ec = []
    for _, row in ec.iterrows():
        all_ec.append(row["ec"])

    counts = dict()
    for i in all_ec:
        counts[i] = dict()

    for _, row in taxons.iterrows():
        raw_ec = EC_RE.search(str(row["readID"]))
        if not raw_ec:
            continue

        read_ec = raw_ec.group(1)
        read_ec = read_ec.split("_")[0]

        tax_id = str(row["taxID"])

        if read_ec not in all_ec:
            continue

        if tax_id not in counts[read_ec]:
            counts[read_ec][tax_id] = 0
        counts[read_ec][tax_id] += 1

    # Create matrix EC (rows) x taxID (columns)
    mat = pd.DataFrame.from_dict(counts, orient="index").fillna(0).astype(int)
    mat.index.name = "EC"

    # Write to Excel
    mat.to_excel(args.output)

if __name__ == "__main__":
    main()

