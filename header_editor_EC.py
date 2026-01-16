#!/usr/bin/env python3
import argparse
import re
from Bio import SeqIO

bracket = re.compile(r"\[([^=\]]+)=([^\]]+)\]")

def fields(desc: str) -> dict:
    return dict(bracket.findall(desc))

def main():
    p = argparse.ArgumentParser(description="Adds EC annotations to FASTA file header header (only keeps annotated results")
    p.add_argument("--gbff", required=True, help="GenBank file directory")
    p.add_argument("--cds", required=True, help="CDS file directory")
    p.add_argument("--out", required=True, help="Output FASTA file directory")
    args = p.parse_args()

    # 1) Build EC lookup from GBFF
    ec_by_pid, ec_by_lt = {}, {}

    for rec in SeqIO.parse(args.gbff, "genbank"):
        for f in rec.features:
            if f.type != "CDS":
                continue
            q = f.qualifiers
            ecs = q.get("EC_number")
            if not ecs:
                continue
            ec = ";".join(ecs)
            pid = (q.get("protein_id") or [""])[0]
            lt = (q.get("locus_tag") or [""])[0]
            if pid and pid not in ec_by_pid:
                ec_by_pid[pid] = ec
            if lt and lt not in ec_by_lt:
                ec_by_lt[lt] = ec

    # 2) Write only sequences with EC
    kept = 0
    total = 0
    with open(args.out, "w") as out:
        for r in SeqIO.parse(args.cds, "fasta"):
            total += 1
            info = fields(r.description)
            pid = info.get("protein_id", "")
            lt = info.get("locus_tag", "")

            ec = ec_by_pid.get(pid) or ec_by_lt.get(lt)
            if not ec:
                continue  # skip no-EC

            base = pid or r.id
            r.id = f"{base}|EC={ec}"  # EC stays in first token
            r.name = r.id
            r.description = ""
            SeqIO.write(r, out, "fasta")
            kept += 1

    print(f"Total CDS: {total}")
    print(f"Kept with EC: {kept}")
    print(f"Wrote: {args.out}")

if __name__ == "__main__":
    main()