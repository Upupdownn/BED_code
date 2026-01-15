#!/usr/bin/env python3
"""
Extract End-Motif (k-mer) frequencies from fragment TSV file.

Input: TSV file with columns: chr start end mapq strand.
Output: TSV file with kmer, count, freq.
"""

import argparse
import logging
import pandas as pd
import py2bit
from pathlib import Path
from multiprocessing import Pool
import math
import gzip


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("frag_tsv_file", help="Input Fragment TSV file")
    parser.add_argument("tb_file", help="Reference genome 2bit file")
    parser.add_argument("out_file", help="Output TSV file for k-mer frequencies")
    parser.add_argument("--len_min", type=int, default=20, help="Minimum fragment length [default: 20]")
    parser.add_argument("--len_max", type=int, default=600, help="Maximum fragment length [default: 600]")
    parser.add_argument("--qmin", type=int, default=30, help="Minimum mapping quality [default: 30]")
    parser.add_argument("--k", type=int, default=4, help="k-mer length [default: 4]")
    parser.add_argument("--no-autosomes", action="store_true", help="Disable autosome-only filtering (default: enabled)")
    parser.add_argument("-p", "--processes", type=int, default=20, help="Number of processes [default: 20]")
    return parser.parse_args()

class Utils():
    @staticmethod
    def kmer_list(k=4) -> list:
        """Generate all possible k-mers recursively."""
        bases = "ACGT"
        if k == 1:
            return list(bases)
        else:
            kmers = []
            for k_minus_mer in Utils.kmer_list(k - 1):
                for base in bases:
                    kmers.append(k_minus_mer + base)
            return kmers
        
    @staticmethod  
    def kmer_counts(k=4):
        """Initialize k-mer count dictionary with zeros."""
        return {kmer: 0 for kmer in Utils.kmer_list(k)}
    
    @staticmethod
    def chr_str2int_map():
        """Mapping for chromosome names (supports both '1' and 'chr1')."""
        return {f"{i}": i for i in range(1, 23)} | {f"chr{i}": i for i in range(1, 23)}
    
    @staticmethod
    def filter_frag_in_df(df: pd.DataFrame, len_min=20, len_max=600, qmin=30, autosomes=True):
        """Filter fragments by length, mapq, and optionally autosomes."""
        df['chr'] = df['chr'].astype(str)
        if autosomes:
            df = df.loc[df['chr'].isin(Utils.chr_str2int_map()), :]

        df = df.loc[df['mapq'] >= qmin, :]
        lens = df['end'] - df['start']
        df = df.loc[(lens >= len_min) & (lens <= len_max)]
        return df
    
    @staticmethod
    def get_seq(tbfile, chrom:str, start, end, close=False, check=False):
        """Safely fetch sequence."""
        target_chrom = chrom
        try:
            if hasattr(tbfile, 'chroms'):
                available = tbfile.chroms()
            elif hasattr(tbfile, 'references'):
                available = tbfile.references
            else:
                available = []

            if available and chrom not in available:
                if chrom.startswith("chr"):
                    alt_chrom = chrom[3:]
                else:
                    alt_chrom = f"chr{chrom}"
                
                if alt_chrom in available:
                    target_chrom = alt_chrom
        except Exception:
            pass

        try:
            seq = tbfile.sequence(target_chrom, start, end)
        except Exception as e:
            seq = "N" * (end - start)

        if check and (len(seq) != end - start or 'N' in seq):
            seq = None
        
        if close: 
            tbfile.close()

        return seq
    
    @staticmethod
    def rvs_cplmt_seq(seq: str) -> str:
        """Return the reverse complement of the sequence."""
        if seq is None: return None
        
        reversed = seq[-1::-1]
        pair_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        res = ''.join(pair_dict[base] for base in reversed)
        return res
    
def count_end_motifs_in_df(df: pd.DataFrame, tbfile_path, k=4) -> dict:
    """Count k-mers at both 5' ends of fragments (forward and reverse complement)."""
    kmer_counts = Utils.kmer_counts(k)
    if df.empty: return kmer_counts

    tbfile = py2bit.open(tbfile_path, 'r')
    for chrom, start, end, _, _ in df.itertuples(index=False):
        pos = start, start + k, end - k, end
        r1_kmer = Utils.get_seq(tbfile, chrom, *pos[:2], close=False, check=True)
        r2_kmer = Utils.rvs_cplmt_seq(Utils.get_seq(tbfile, chrom, *pos[2:], close=False, check=True))
        if r1_kmer is None or r2_kmer is None: continue
        kmer_counts[r1_kmer] += 1
        kmer_counts[r2_kmer] += 1
    tbfile.close()
    return kmer_counts

def count_chunk_from_tsv(args):
    """
    Multi-processing entry function:
    Read a specific line range from TSV file, filter, then count k-mers.
    """
    tsv_file, start_line, end_line, tbfile_path, k, len_min, len_max, qmin, autosomes = args

    # Read specified line range (skip header if start_line==0)
    skip_rows = 1 if start_line == 0 else start_line + 1                # +1 to skip header
    nrows = end_line - start_line if end_line is not None else None

    df_chunk = pd.read_table(
        tsv_file,
        header=None if start_line > 0 else 0,
        skiprows=skip_rows if start_line > 0 else 0,
        nrows=nrows,
        names=['chr', 'start', 'end', 'mapq', 'strand'],
        dtype={'chr': str}
    )
    logger.info(f"Processing lines {start_line+1} to {end_line if end_line else 'end'} ({len(df_chunk)} rows)")

    df_filtered = Utils.filter_frag_in_df(df_chunk, len_min, len_max, qmin, autosomes)

    return count_end_motifs_in_df(df_filtered, tbfile_path, k=k)
    
def extract_edm_features(tsv_file, tb_file, out_file, len_min=20, len_max=600, qmin=30, autosomes=True, k=4, processes=20):
    """Main function: Load, split, count k-mers in parallel, compute frequencies, save."""
    is_gz = str(tsv_file).lower().endswith('.gz')
    opener = gzip.open if is_gz else open

    # Count total lines (excluding header)
    logger.info(f"Loading TSV file to get total lines: {tsv_file}")
    total_lines = 0
    with opener(tsv_file, 'rt') as f:
        total_lines = sum(1 for _ in f) - 1
    logger.info(f"Total data rows (excluding header): {total_lines:,}")

    if total_lines <= 0:
        logger.warning("No data rows found. Exit.")
        return
    
    # Split into chunks of ~0.5 million rows each
    chunk_size = 500000
    num_chunks = math.ceil(total_lines / chunk_size)
    tasks = []
    for i in range(num_chunks):
        start = i * chunk_size
        end = min(start + chunk_size, total_lines)
        tasks.append((tsv_file, start, end, tb_file, k, len_min, len_max, qmin, autosomes))
    logger.info(f"Splitting into {num_chunks} chunks for {processes} processes")

    # Parallel counting
    with Pool(processes=processes) as pool:
        results = pool.map(count_chunk_from_tsv, tasks)

    # Merge all counts
    logger.info("Merging counts from all processes...")
    total_counts = Utils.kmer_counts(k)
    for res in results:
        for kmer, cnt in res.items():
            total_counts[kmer] += cnt

    # Compute frequencies and save
    logger.info("Calculating frequencies...")
    result_df = pd.DataFrame.from_dict(total_counts, orient='index', columns=['count'])
    total = result_df['count'].sum()
    result_df['freq'] = result_df['count'] / total if total > 0 else 0
    result_df.index.name = 'kmer'

    logger.info(f"Saving results to: {out_file}")
    result_df.to_csv(out_file, sep='\t', header=True, index=True)
    logger.info(f"Completed. Total fragments processed: {total_lines:,}, Total counts: {total:,}")

def main():
    args = parse_args()

    extract_edm_features(
        tsv_file=args.frag_tsv_file,
        tb_file=args.tb_file,
        out_file=args.out_file,
        len_min=args.len_min,
        len_max=args.len_max,
        qmin=args.qmin,
        autosomes=not args.no_autosomes,  # Default is True (enabled)
        k=args.k,
        processes=args.processes
    )

if __name__ == "__main__":
    main()

