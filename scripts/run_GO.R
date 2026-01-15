#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# GO Enrichment Analysis for Human Genes (SYMBOL → ENTREZID → enrichGO)
#
# This script performs Gene Ontology (GO) enrichment analysis using clusterProfiler
# on a list of gene symbols provided in a tab-separated input file.
#
# Input:  TSV file with at least a "Gene" column containing HGNC symbols
# Output: TSV file with enriched GO terms (BP by default)
#
# Usage:
#   Rscript run_GO.R input_genes.tsv output_go.tsv
#   Rscript run_GO.R input_genes.tsv output_go.tsv --ont MF --pvalue 0.01
#
# Dependencies:
#   clusterProfiler, org.Hs.eg.db, argparser
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(argparser)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# Parse command-line arguments
p=arg_parser("Perform GO enrichment analysis on a list of human gene symbols.")
p=add_argument(p,"input_file",help="Input TSV file containing gene symbols")
p=add_argument(p,"output_file",help="Output TSV file for enriched GO terms")
p=add_argument(p,"--gene_col",default="Gene",help="Name of the gene symbol column")
p=add_argument(p,"--ont",default="BP",help="GO ontology: BP, MF, CC, or ALL")
p=add_argument(p,"--pvalue",type="numeric",default=0.05,help="p-value cutoff")
p=add_argument(p,"--qvalue",type="numeric",default=0.05,help="q-value(FDR) cutoff")
p=add_argument(p,"--pAdjustMethod",default="BH",help="p-value adjustment method")
args=parse_args(p)

# Main function
main <- function() {
  cat("GO Enrichment Analysis\n")
  cat("======================\n")
  cat("Input file:  ", args$input_file, "\n")
  cat("Gene column: ", args$gene_col, "\n")
  cat("Output file: ", args$output_file, "\n")
  cat("Ontology:    ", args$ont, "\n")
  cat("p-value cut: ", args$pvalue, "\n")
  cat("q-value cut: ", args$qvalue, "\n\n")

  # Read input genes
  if (!file.exists(args$input_file)) 
  {
    stop("Error: Input file not found: ", args$input_file)
  }

  df <- read.table(args$input_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

  if (!args$gene_col %in% colnames(df)) 
  {
    stop("Error: Input file must contain a gene column named")
  }
  colnames(df)[colnames(df) == args$gene_col] <- "Gene"
  
  genes <- unique(df$Gene)
  genes <- genes[genes != "" & !is.na(genes)]
  cat("Number of unique input genes:", length(genes), "\n")

  if (length(genes) == 0) 
  {
    stop("Error: No valid genes found in input file.")
  }

  # SYMBOL to ENTREZID conversion
  cat("Converting gene symbols to ENTREZ IDs...\n")
  gene_df <- bitr(genes,
                  fromType   = "SYMBOL",
                  toType     = "ENTREZID",
                  OrgDb      = org.Hs.eg.db,
                  drop       = TRUE)

  if (nrow(gene_df) == 0) 
  {
    stop("Error: No genes could be mapped to ENTREZ IDs.")
  }

  entrez_ids <- unique(gene_df$ENTREZID)
  cat("Number of successfully mapped ENTREZ IDs:", length(entrez_ids), "\n\n")

  # Perform GO enrichment
  cat("Running enrichGO...\n")
  ego <- enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = args$ont,
    pAdjustMethod = args$pAdjustMethod,
    pvalueCutoff  = args$pvalue,
    qvalueCutoff  = args$qvalue,
    readable      = TRUE,
  )

  # Save results
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) 
  {
    cat("Warning: No significant GO terms found.\n")
    write.table(data.frame(), file=args$output_file, sep= "\t", quote= FALSE, row.names=FALSE)
  } 
  else 
  {
    cat("Number of enriched GO terms:", nrow(as.data.frame(ego)), "\n")
    write.table(as.data.frame(ego),
                file      = args$output_file,
                sep       = "\t",
                quote     = FALSE,
                row.names = FALSE)
    cat("Results saved to:", args$output_file, "\n")
  }
}

# Run
main()