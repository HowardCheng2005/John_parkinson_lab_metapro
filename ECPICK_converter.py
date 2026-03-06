# run_ecpick.py
import argparse
from ecpick import ECPICK

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--fasta", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    args = ap.parse_args()

    model = ECPICK()
    # Writes args.outdir/result.csv (or result.json if you choose)
    model.predict_fasta(fasta_path=args.fasta, output_path=args.outdir)

if __name__ == "__main__":
    main()