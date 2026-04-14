#!/usr/bin/env python3
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq

bracket = re.compile(r"\[([^=\]]+)=([^\]]+)\]")

def fields(desc: str) -> dict:
    return dict(bracket.findall(desc))

def main():
    p = argparse.ArgumentParser(description="Adds EC annotations to FASTA file header header (only keeps annotated results")
    p.add_argument("--gbff", required=True, help="GenBank file directory")
    p.add_argument("--cds", required=True, help="CDS file directory")
    p.add_argument("--nucleotide_out", required=True, help="Output FASTA file directory for nucleotide annotations")
    p.add_argument("--translation_out", required=True, help="Output FASTA file directory for translation")
    args = p.parse_args()

    # 1) Build EC lookup from GBFF
    ec_by_pid, ec_by_lt = {}, {}
    translation_by_pid, translation_by_lt = {}, {}

    for rec in SeqIO.parse(args.gbff, "genbank"):
        for f in rec.features:
            if f.type != "CDS":
                continue
            q = f.qualifiers
            translation = (q.get("translation") or [""])[0]
            pid = (q.get("protein_id") or [""])[0]
            lt = (q.get("locus_tag") or [""])[0]

            if translation:
                if pid and pid not in translation_by_pid:
                    translation_by_pid[pid] = translation
                if lt and lt not in translation_by_lt:
                    translation_by_lt[lt] = translation

            ecs = q.get("EC_number")
            if not ecs:
                continue
            ec = ";".join(ecs)

            if pid and pid not in ec_by_pid:
                ec_by_pid[pid] = ec
            if lt and lt not in ec_by_lt:
                ec_by_lt[lt] = ec

    # 2) Write sequences with EC number (-1 if otherwise)
    nucleotide_with_ec = 0
    nucleotide_total = 0
    with open(args.nucleotide_out, "w") as out:
        for i, r in enumerate(SeqIO.parse(args.cds, "fasta")):
            nucleotide_total += 1
            info = fields(r.description)
            pid = info.get("protein_id", "")
            lt = info.get("locus_tag", "")

            ec = ec_by_pid.get(pid) or ec_by_lt.get(lt) or "-1"
            if ec != "-1":
                nucleotide_with_ec += 1  # if EC

            base = pid or r.id
            r.id = f"{i:09d}|{base}|EC={ec}"  # EC stays in first token
            r.name = r.id
            r.description = ""
            SeqIO.write(r, out, "fasta")

    translation_with_ec = 0
    translation_total = 0

    with open(args.translation_out, "w") as out:
        for j, r in enumerate(SeqIO.parse(args.cds, "fasta")):
            translation_total += 1
            info = fields(r.description)
            pid = info.get("protein_id", "")
            lt = info.get("locus_tag", "")

            ec = ec_by_pid.get(pid) or ec_by_lt.get(lt) or "-1"
            if ec != "-1":
                translation_with_ec += 1
            translation = translation_by_pid.get(pid) or translation_by_lt.get(lt)
            if not translation:
                continue

            base = pid or r.id
            r.seq = Seq(translation)
            r.id = f"{j:09d}|{base}|EC={ec}"
            r.name = r.id
            r.description = ""
            SeqIO.write(r, out, "fasta")

    print(f"Total CDS: {nucleotide_total}")
    print(f"With EC: {nucleotide_with_ec}")
    print(f"Wrote: {args.nucleotide_out}")

    print(f"Total CDS: {translation_total}")
    print(f"With EC: {translation_with_ec}")
    print(f"Wrote: {args.translation_out}")

if __name__ == "__main__":
    main()
