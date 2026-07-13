############################################################
# Project: TCGA-RNAseq-Analysis
# Script : 05_gene_annotation.R
# Author : Zeinab Rohani
#
# Description:
# Annotate differentially expressed genes by mapping
# Ensembl gene IDs to official HGNC gene symbols.
############################################################

##==========================================================
## Set Working Directory
##==========================================================

setwd("D:/TCGA/")

##==========================================================
## Load Gene Annotation
##==========================================================

gene_annotation <- read.delim(
  "results/gene_annotation.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

##==========================================================
## Load Differential Expression Results
##==========================================================

deg_results <- read.delim(
  "results/differential_expression_results.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

##==========================================================
## Prepare Annotation Table
##==========================================================

rownames(gene_annotation) <- gene_annotation$ensembl_gene_id

##==========================================================
## Remove Version Numbers (Optional)
##==========================================================

deg_results$Ensembl_ID <- sub(
  "\\..*$",
  "",
  deg_results$Ensembl_ID
)

##==========================================================
## Map Gene Symbols
##==========================================================

deg_results$Gene_Symbol <- gene_annotation[
  deg_results$Ensembl_ID,
  "external_gene_name"
]

##==========================================================
## Reorder Columns
##==========================================================

deg_results <- deg_results[
  ,
  c(
    "Ensembl_ID",
    "Gene_Symbol",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
]

##==========================================================
## Export Annotated Results
##==========================================================

write.table(
  deg_results,
  file = "results/annotated_DEGs.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

##==========================================================
## Export Significant Annotated DEGs
##==========================================================

significant_deg <- subset(
  deg_results,
  abs(logFC) >= 1 &
    adj.P.Val < 0.05
)

write.table(
  significant_deg,
  file = "results/significant_annotated_DEGs.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

##==========================================================
## Summary
##==========================================================

cat("\n")
cat("=========================================\n")
cat("Gene Annotation Completed\n")
cat("-----------------------------------------\n")
cat("Annotated genes :", nrow(deg_results), "\n")
cat("Significant DEGs:", nrow(significant_deg), "\n")
cat("-----------------------------------------\n")
cat("Results exported to results/\n")
cat("=========================================\n")
