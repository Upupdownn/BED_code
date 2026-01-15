#!/usr/bin/env python3
"""
BAM File to Fragment TSV File Converter.
TSV File Columns: chr | start | end | mapq | strand (0-based start, 1-based end)
"""


import argparse
import pysam
import pandas as pd
import logging
from multiprocessing import Pool
from pathlib import Path


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("out_tsv_file", help="Output TSV file")
    parser.add_argument("-p", "--processes", type=int, default=20, help="Number of processes [default: 20]")
    parser.add_argument("--mapq-min", type=int, default=30, help="Minimum mapping quality [default: 30]")

    return parser.parse_args()

def is_invalid_read(read: pysam.AlignedSegment, mapq_min: int) -> bool:
    if ((not read.is_paired)
        or read.is_secondary
        or read.is_unmapped
        or read.mate_is_unmapped
        or read.mapping_quality < mapq_min
        or read.is_qcfail
        or read.is_duplicate
        or read.is_supplementary
        or (not read.is_proper_pair)
        or (read.is_reverse and read.mate_is_reverse)
        ):
        return True
    else:
        return False
    
def process_chrom(args):
    bam_path, chrom, mapq_min = args

    rows = []
    qname_mapq_dict = dict()
    bad_count = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam_file:
        for read in bam_file.fetch(chrom):
            if is_invalid_read(read, mapq_min):
                bad_count += 1
                continue
            
            if read.query_name not in qname_mapq_dict:
                qname_mapq_dict[read.query_name] = read.mapping_quality
                continue

            ref_name = read.reference_name
            seq_len = abs(read.template_length)
            if read.is_read1:
                mapped_strand = '-' if read.is_reverse else '+'         
            else:
                mapped_strand = '+' if read.is_reverse else '-'

            mapped_quality = min([read.mapping_quality, qname_mapq_dict.pop(read.query_name, None)])

            if read.is_reverse:
                ref_end = read.reference_end
                ref_start = ref_end - seq_len
            else:
                ref_start = read.reference_start
                ref_end = ref_start + seq_len

            rows.append([ref_name, ref_start, ref_end, mapped_quality, mapped_strand])

        bad_count += len(qname_mapq_dict)

    df = pd.DataFrame(rows, columns=['chr', 'start', 'end', 'mapq', 'strand'])
    df.sort_values(by=['start', 'end'], inplace=True)

    return df, bad_count

def main():
    args = parse_args()

    bam_path = Path(args.bam_file)
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if not bam.has_index():
            logger.error(f"Error: BAM index (.bai) missing for {bam_path}")
        chroms = list(bam.header.references)

    logger.info(f"Processing {len(chroms)} chromosomes with {args.processes} processes")

    arg_list = [(str(bam_path), chrom, args.mapq_min) for chrom in chroms]
    df_list = []
    total_bad = 0

    with Pool(args.processes) as pool:
        for df, bad in pool.imap(process_chrom, arg_list):
            df_list.append(df)
            total_bad += bad

    logger.info(f"Total bad reads filtered: {total_bad:,}")

    final_df = pd.concat(df_list, ignore_index=True)
    final_df.to_csv(args.out_tsv_file, sep='\t', header=True, index=False)
    logger.info(f"Output saved to {args.out_tsv_file} ({len(final_df):,} fragments)")


if __name__ == "__main__":
    main()
