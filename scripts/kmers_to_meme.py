#!/home/shengxinwei/miniconda3/envs/finale/bin/python3.10
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("kmer_list_file")
parser.add_argument("meme_out_file")
cmd_args = parser.parse_args()

kmers_file = cmd_args.kmer_list_file   # 每行一个6-mer, 例如 "ACGTGA"
out_meme = cmd_args.meme_out_file

alphabet = ["A","C","G","T"]

def kmer_to_pwm(kmer):
    '''给定一个kmer，输出PWM矩阵'''
    rows = []
    for ch in kmer:
        row = [ "1.0" if ch == b else "0.0" for b in alphabet ]
        rows.append(row)
    return rows

with open(kmers_file) as fin, open(out_meme,"w") as fo:
    fo.write("MEME version 4\n\nALPHABET= ACGT\n\n")
    fo.write("strands: + -\n\n")
    fo.write("Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
    for line in fin:
        k = line.strip().upper()
        if not k: continue
        fo.write(f"MOTIF {k}\n")
        # alength:字母表长度，w:motif长度，nsites:motif由多少序列组成，E:motif期望值(均为注释)
        fo.write(f"letter-probability matrix: alength= 4 w= {len(k)} nsites= 1 E= 0\n")
        pwm = kmer_to_pwm(k)
        for row in pwm:
            fo.write(" ".join(row) + "\n")
        fo.write("\n")

print("Wrote", out_meme)
