#!/usr/bin/env python3
"""
Map TOMTOM output (Query -> Target) to TRRUST2 TF-Gene regulations.

This script:
1. Reads TOMTOM results (k-mer -> TF name)
2. Maps each TF to its regulated genes from TRRUST2
3. Outputs a tidy table: kmer | TF | Gene | Regulation | PMID
4. Prints a summary of unique k-mers, TFs, and genes

Usage:
    python tomtom_to_trrust.py tomtom.tsv trrust2.tsv output.tsv
"""

import argparse
from pathlib import Path
import pandas as pd


def load_trrust2(filepath: Path) -> pd.DataFrame:
    """Load TRRUST2 data with proper column names."""
    df = pd.read_table(filepath, header=None, usecols=[0, 1, 2, 3], names=["TF", "Gene", "Regulation", "PMID"])

    return df

def load_tomtom(tomtom_path: Path) -> pd.DataFrame:
    """Load TOMTOM output, keeping only Query_ID (k-mer) and Target_ID (TF)."""
    df = pd.read_table(tomtom_path, comment="#", usecols=["Query_ID", "Target_ID"], skip_blank_lines=True)

    # Clean TF name: remove everything after the first dot (e.g., "TF.1" -> "TF")
    df["Target_ID"] = df["Target_ID"].str.split(".").str[0].str.strip()
    df = df.rename(columns={"Query_ID": "kmer", "Target_ID": "TF"})
    return df

def map_to_trrust(tomtom_df: pd.DataFrame, trrust_df: pd.DataFrame) -> pd.DataFrame:
    """
    Join TOMTOM k-mers with TRRUST2 regulations.
    Returns a tidy DataFrame with kmer repeated for each TF-Gene pair.
    """
    # Merge on TF name (inner join: only keep matched TFs)
    merged = tomtom_df.merge(trrust_df, on="TF", how="inner")

    # Reorder columns for clarity
    merged = merged[["kmer", "TF", "Gene", "Regulation", "PMID"]]

    # Sort by kmer for consistent output
    merged = merged.sort_values("kmer")

    return merged


def print_summary(df: pd.DataFrame):
    """Print summary statistics: unique k-mers, TFs, genes per k-mer and overall."""
    print("\nSummary:")
    print("kmer\t#TF\t#Gene")

    for kmer in sorted(df["kmer"].unique()):
        sub = df[df["kmer"] == kmer]
        n_tf = sub["TF"].nunique()
        n_gene = sub["Gene"].nunique()
        print(f"{kmer}\t{n_tf}\t{n_gene}")

    total_kmer = df["kmer"].nunique()
    total_tf = df["TF"].nunique()
    total_gene = df["Gene"].nunique()
    print(f"\nTotal: motif={total_kmer}\tTF={total_tf}\tGene={total_gene}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("tomtom_file", type=Path, help="TOMTOM output file (tsv)")
    parser.add_argument("trrust2_file", type=Path, help="TRRUST2 raw data file (human.tsv)")
    parser.add_argument("output_file", type=Path, help="Output tidy table (tsv)")

    args = parser.parse_args()

    # Validate input files exist
    if not args.tomtom_file.is_file():
        parser.error(f"TOMTOM file not found: {args.tomtom_file}")
    if not args.trrust2_file.is_file():
        parser.error(f"TRRUST2 file not found: {args.trrust2_file}")

    print(f"Loading TOMTOM results from: {args.tomtom_file}")
    tomtom_df = load_tomtom(args.tomtom_file)

    print(f"Loading TRRUST2 from: {args.trrust2_file}")
    trrust_df = load_trrust2(args.trrust2_file)

    print("Mapping k-mers to regulations...")
    result_df = map_to_trrust(tomtom_df, trrust_df)

    if result_df.empty:
        print("Warning: No matches found between TOMTOM and TRRUST2.")
    else:
        print(f"Found {len(result_df)} regulation records.")
        result_df.to_csv(args.output_file, sep="\t", index=False)
        print(f"Saved tidy table to: {args.output_file}")

    print_summary(result_df)


if __name__ == "__main__":
    main()