#!/home/shengxinwei/miniconda3/envs/finale/bin/python3.10

import argparse
import pandas as pd

trrust2_fpath = "/home/shengxinwei/biodata/Files/TRRUST2/trrust_rawdata.human.tsv"  # TF和gene对应

parser = argparse.ArgumentParser()
parser.add_argument("tomtom_out_file")
parser.add_argument("tidy_out_file")
cmd_args = parser.parse_args()

def summary(df: pd.DataFrame):
    print(f"kmer\tn_TF\tn_Gene")
    for kmer in df['kmer'].unique():
        sub_df = df.loc[df["kmer"] == kmer, :]
        n_TF = len(sub_df["TF"].unique())
        n_Gene = len(sub_df["Gene"].unique())
        print(f"{kmer}\t{n_TF}\t{n_Gene}")
    n_kmer = len(df["kmer"].unique())
    n_TF = len(df["TF"].unique())
    n_Gene = len(df["Gene"].unique())
    print(f"motif:{n_kmer}\tTF:{n_TF}\tGene:{n_Gene}")
    pass

tomtom_file = cmd_args.tomtom_out_file
out_file = cmd_args.tidy_out_file
tomtom_df = pd.read_table(tomtom_file, usecols=["Query_ID", "Target_ID"], comment='#', skip_blank_lines=True)
trrust2_df = pd.read_table(trrust2_fpath, header=None, index_col=None)
trrust2_df.columns = ["TF", "Gene", "Regulation", "PMID"]

result_df = pd.DataFrame()
for kmer, TF_name in tomtom_df.itertuples(index=False):
    TF_name: str = TF_name.strip().split('.')[0]
    sub_df = trrust2_df.loc[trrust2_df["TF"] == TF_name, :].iloc[:, :3].copy()
    if (len(sub_df) == 0): continue
    sub_df.insert(0, 'kmer', kmer)
    result_df = pd.concat([result_df, sub_df], axis=0)

result_df.sort_values(by=["kmer"], inplace=True)
result_df.to_csv(out_file, sep='\t', header=True, index=False)
summary(result_df)


