#!/usr/bin/env python3
"""
Compute reference background and pilot subset from EDM feature matrix.

Input:
  - edm_feature_file: TSV with samples (rows) x 256 motifs (columns)
  - info_file: TSV with at least sample id and label columns

Output (optional):
  - --ref_out: npy file of mean frequency for label=0 samples (reference background)
  - --pilot_out: TSV file of pilot subset (20% stratified by label)
"""

import argparse
import logging
import pandas as pd
import numpy as np
from pathlib import Path
import os
from sklearn.model_selection import train_test_split

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("edm_feature_file", help="TSV file: samples x 256 motif frequencies")
    parser.add_argument("info_file", help="TSV file with at least id and label columns")
    parser.add_argument("--id_col", default="id", help="Column name for sample ID [default: id]")
    parser.add_argument("--label_col", default="label", help="Column name for label (0/1) [default: label]")
    parser.add_argument("--ref_out", default=None, help="Output path for reference background (.npy)")
    parser.add_argument("--pilot_out", default=None, help="Output path for pilot subset features (.tsv)")
    parser.add_argument("--pilot_ratio", type=float, default=0.2, help="Pilot subset ratio [default: 0.2]")
    parser.add_argument("--random_state", type=int, default=42, help="Random seed [default: 42]")
    return parser.parse_args()


def main():
    args = parse_args()

    # Load feature matrix (samples x motifs)
    logger.info(f"Loading feature matrix: {args.edm_feature_file}")
    features = pd.read_table(args.edm_feature_file, index_col=0)

    # Load metadata (id and label)
    logger.info(f"Loading metadata: {args.info_file}")
    metadata = pd.read_table(args.info_file, index_col=None)

    # Ensure sample IDs match
    common_ids = features.index.intersection(metadata[args.id_col])
    if len(common_ids) == 0:
        raise ValueError("No common sample IDs between features and metadata")
    
    features = features.loc[common_ids]
    metadata = metadata.set_index(args.id_col).loc[common_ids]

    logger.info(f"Get {len(common_ids)} samples")

    # Task 1: Compute reference background (mean frequency of label=0 samples)
    if args.ref_out:
        control_ids = metadata[metadata[args.label_col] == 0].index
        if len(control_ids) == 0:
            raise ValueError("No control samples (label=0) found")
        
        ref_bg = features.loc[control_ids].mean(axis=0).values
        os.makedirs(os.path.dirname(args.ref_out), exist_ok=True)
        np.save(args.ref_out, ref_bg)
        logger.info(f"Reference background (label=0 mean) saved to: {args.ref_out} (shape: {ref_bg.shape})")

    # Task 2: Stratified sampling for pilot set (20% by default)
    if args.pilot_out:
        # Stratified split on features and labels
        pilot_df, _ = train_test_split(
            features,
            test_size=1 - args.pilot_ratio,
            stratify=metadata[args.label_col],
            random_state=args.random_state
        )
        os.makedirs(os.path.dirname(args.pilot_out), exist_ok=True)
        pilot_df.to_csv(args.pilot_out, sep='\t')
        logger.info(f"Pilot subset ({len(pilot_df)} samples, {args.pilot_ratio:.1%}) saved to: {args.pilot_out}")
    logger.info("Completed successfully.")


if __name__ == "__main__":
    main()