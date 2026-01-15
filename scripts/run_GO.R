#!/home/shengxinwei/miniconda3/envs/finale/bin/Rscript

suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))

#  解析命令行参数
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: run_GO.R <input_file> <output_file>\n")
}

input_file  <- args[1]
output_file <- args[2]

cat("Input file: ", input_file, "\n")
cat("Output file:", output_file, "\n")

# 读取输入文件(该文件需要有Gene列)
df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- unique(df[["Gene"]])

cat("Number of input genes:", length(genes), "\n")

# 基因 ID 转换 SYMBOL → ENTREZ
gene_df <- bitr(
    genes,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
)

if (nrow(gene_df) == 0) {
    stop("ERROR: No valid gene IDs found after converting SYMBOL → ENTREZID.\n")
}

entrez_ids <- unique(gene_df$ENTREZID)
cat("Number of ENTREZ IDs:", length(entrez_ids), "\n")

# GO 富集分析（BP）
# GO注释类别，BP(生物学过程), MF(分子功能), CC(细胞组成)
ego <- enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",   
    pAdjustMethod = "BH",   # FDR矫正
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
)

cat("GO terms found:", nrow(as.data.frame(ego)), "\n")

# 保存结果
write.table(as.data.frame(ego), file=output_file, sep="\t", quote=FALSE, row.names=FALSE)

cat("GO enrichment results saved to:", output_file, "\n")
