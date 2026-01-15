#!/usr/bin/env python3
"""
Convert a list of k-mers (one per line) into a MEME motif file format.
Each k-mer is treated as a separate motif with a deterministic PWM (1.0 at the matching base, 0.0 elsewhere).

Usage:
    ./kmer_to_meme.py kmers.txt motifs.meme
"""

import argparse
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("kmer_list_file", type=Path, help="Input file: one uppercase k-mer per line")
    parser.add_argument("meme_out_file", type=Path, help="Output MEME file path")

    return parser.parse_args()

def kmer_to_pwm(kmer: str, alphabet: str = "ACGT") -> list[list[str]]:
    """Convert a single k-mer into a PWM matrix."""
    return [["1.0" if base == ch else "0.0" for base in alphabet] for ch in kmer]


def main():
    args = parse_args()

    alphabet = "ACGT"
    header = [
        "MEME version 4",
        "",
        "ALPHABET= ACGT",
        "",
        "strands: + -",
        "",
        "Background letter frequencies",
        "A 0.25 C 0.25 G 0.25 T 0.25",
        "",
    ]

    with args.kmer_list_file.open() as fin, args.meme_out_file.open("w") as fo:
        # Write header once
        fo.write("\n".join(header) + "\n")

        motif_count = 0
        for line in fin:
            kmer = line.strip().upper()
            if not kmer or len(kmer) == 0:
                continue

            if not all(base in alphabet for base in kmer):
                logging.info(f"Warning: skipping invalid k-mer: {kmer}")
                continue

            motif_count += 1

            # Write motif block
            fo.write(f"MOTIF {kmer}\n")
            fo.write(f"letter-probability matrix: alength= 4 w= {len(kmer)} nsites= 1 E= 0\n")

            pwm = kmer_to_pwm(kmer, alphabet)
            for row in pwm:
                fo.write(" ".join(row) + "\n")
            fo.write("\n")

    logging.info(f"Done. Wrote {motif_count} motifs to {args.meme_out_file}")


if __name__ == "__main__":
    main()