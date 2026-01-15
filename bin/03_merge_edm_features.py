#!/usr/bin/env python3
"""
Merge individual sample EDM feature files into a single matrix.

Input: Directory containing one TSV file per sample (k-mer frequencies)
Output:
  - Merged feature matrix (samples x k-mers)
  - Optional: Motif Diversity Score (MDS) for each sample
"""

import argparse
import logging
import pandas as pd
import numpy as np
from pathlib import Path
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("input_dir", help="Directory containing per-sample End-Motif Frequency TSV files")
    parser.add_argument("--merged_output", default=None, help="Path to save merged feature matrix (optional)")
    parser.add_argument("--mds_output", default=None, help="Path to save Motif Diversity Score (MDS) results (optional)")
    return parser.parse_args()

def compute_mds(x: np.ndarray):
    """Calculate Motif Diversity Score (MDS) as normalized Shannon entropy."""
    x = x.astype(float)
    if x.ndim == 1: x = x.reshape((1, -1))
    res = np.sum(-x * np.log(x, out=np.zeros_like(x), where=(x > 1e-9)), axis=1)
    res /= np.log(x.shape[1])
    if res.shape[0] == 1: res = res[0]
    return res

def merge_edm_features(input_dir: str, merged_output: str = None, mds_output: str = None):
    """Merge all sample EDM feature files into one matrix."""
    input_path = Path(input_dir)
    if not input_path.is_dir():
        raise ValueError(f"Input path is not a directory: {input_dir}")

    # Find all TSV files in the directory
    tsv_files = list(input_path.glob("*.tsv"))
    if not tsv_files:
        raise ValueError(f"No TSV files found in directory: {input_dir}")

    logger.info(f"Found {len(tsv_files)} sample feature files")

    # Read and merge
    merged_df = None
    sample_ids = []

    for file_path in tsv_files:
        sample_id = file_path.stem  # Use filename without extension as sample ID

        df = pd.read_table(file_path, index_col=0)  # Assume index = k-mer, column = freq

        # Assume first frequency column (or the only one)
        if 'freq' in df.columns:
            freq_series = df['freq']
        else:
            freq_series = df.iloc[:, 0]  # Take first column as frequency

        freq_series.name = sample_id
        freq_df = pd.DataFrame(freq_series).T  # Make it row vector

        if merged_df is None:
            merged_df = freq_df
        else:
            merged_df = pd.concat([merged_df, freq_df], axis=0)

        sample_ids.append(sample_id)

    # Set index to sample IDs
    merged_df.index = sample_ids
    merged_df.index.name = 'id'

    # Save merged features if requested
    if merged_output:
        Path(merged_output).parent.mkdir(parents=True, exist_ok=True)
        merged_df.to_csv(merged_output, sep='\t')
        logger.info(f"Merged feature matrix saved to: {merged_output}")

    if mds_output is not None:
        mds_scores = compute_mds(merged_df.to_numpy())
        mds_df = pd.DataFrame({'id': sample_ids, 'MDS': mds_scores})
        Path(mds_output).parent.mkdir(parents=True, exist_ok=True)
        mds_df.to_csv(mds_output, sep='\t', header=True, index=False)
        logger.info(f"MDS scores saved to: {mds_output}")

    logger.info(f"Completed. Total samples merged: {len(sample_ids)}")

    return merged_df

def main():
    args = parse_args()

    merge_edm_features(
        input_dir=args.input_dir,
        merged_output=args.merged_output,
        mds_output=args.mds_output
    )

if __name__ == "__main__":
    main()