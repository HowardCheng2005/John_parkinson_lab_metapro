import argparse
import subprocess
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset_root", help="Root folder of dataset")
    parser.add_argument("out_dir", help="Output folder")
    parser.add_argument("header_py", help="path to header_py file")
    args = parser.parse_args()

    dataset_root = Path(args.dataset_root) / "ncbi_dataset" / "data"
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    header_py = Path(args.header_py)

    for dir_temp in dataset_root.iterdir():
        if not dir_temp.is_dir():
            continue

        gbff = next(dir_temp.rglob("*genomic.gbff"), None)
        cds = next(dir_temp.rglob("*cds_from_genomic.fna"), None)

        if gbff is None or cds is None:
            print("missing genomic.gbff or cds_from_genomic.fna file")
            continue

        out_fa = out_dir / f"{dir_temp.name}.fa"

        subprocess.run(
            [sys.executable, args.header_py, "--gbff", str(gbff), "--cds", str(cds), "--out", str(out_fa)],
            check=True
        )

        print ("done: ", dir_temp.name, " at ", out_fa)

    print ("finished!")

if __name__ == "__main__":
    main()