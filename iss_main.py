#!/usr/bin/env python3
import argparse
import configparser
import os
import subprocess
import sys

def install_insilicoseq():
    subprocess.run(
        [sys.executable, "-m", "pip", "install", "--quiet", "insilicoseq"],
        check=True,
    )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input FASTA")
    parser.add_argument("-c", "--config", required=True, help="INI config file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    args = parser.parse_args()

    cfg = configparser.ConfigParser()
    cfg.read(args.config)
    sim = cfg["Simulation"]

    sample_name = sim.get("sample_name", "sample")
    model = sim.get("model", "miseq")
    n_reads = sim.get("n_reads", "1000000")
    abundance = sim.get("abundance", "lognormal")
    sequence_type = sim.get("sequence_type", "metagenomics")
    gc_bias = sim.getboolean("gc_bias", False)
    cpus = sim.get("cpus", "4")
    seed = sim.get("seed", "42")

    os.makedirs(args.output, exist_ok=True)
    out_prefix = os.path.join(args.output, sample_name)

    install_insilicoseq()

    cmd = [
        "iss", "generate",
        "--genomes", args.input,
        "--model", model,
        "--n_reads", n_reads,
        "--abundance", abundance,
        "--sequence_type", sequence_type,
        "--cpus", cpus,
        "--seed", seed,
        "--output", out_prefix,
        "--compress",
    ]
    if gc_bias:
        cmd.append("--gc_bias")

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()