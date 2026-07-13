############################################################
# Project: TCGA-RNAseq-Analysis
# Script : 03_normalization.R
# Author : Zeinab Rohani
#
# Description:
# Filter lowly expressed genes and normalize TCGA RNA-seq
# count data using edgeR and limma-voom.
############################################################

##==========================================================
## Load Required Packages
##==========================================================

library(edgeR)
library(limma)

##==========================================================
## Set Working Directory
##==========================================================

setwd("D:/TCGA/")

##==========================================================
## Load Raw Count Matrix
##==========================================================

counts <- read.delim(
  "results/raw_counts.txt",
  row.names = 1,
  check.names = FALSE
)

##==========================================================
## Load Sample Barcodes
##==========================================================

barcode <- read.delim(
  "barcod/barcod.txt",
  stringsAsFactors = FALSE
)

colnames(counts) <- barcode$barcod

##==========================================================
## Define Experimental Groups
##==========================================================

group <- factor(
  c(
    rep("Normal", 41),
    rep("Cancer", 480)
  )
)

##==========================================================
## Design Matrix
##==========================================================

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

##==========================================================
## Create DGEList Object
##==========================================================

dge <- DGEList(counts = counts)

##==========================================================
## Filter Lowly Expressed Genes
##==========================================================

keep <- filterByExpr(
  dge,
  design = design
)

dge <- dge[keep, , keep.lib.sizes = FALSE]

cat("Genes retained after filtering:", nrow(dge), "\n")

##==========================================================
## TMM Normalization
##==========================================================

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

##==========================================================
## Voom Transformation
##==========================================================

voom_data <- voom(
  dge,
  design,
  plot = TRUE
)

##==========================================================
## Extract Normalized Expression Matrix
##==========================================================

normalized_matrix <- voom_data$E

##==========================================================
## Log2 Transformation
##==========================================================

normalized_matrix <- 2^normalized_matrix
normalized_matrix <- log2(normalized_matrix + 1)

##==========================================================
## Export Normalized Matrix
##==========================================================

write.table(
  normalized_matrix,
  file = "results/normalized_expression_matrix.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

##==========================================================
## Summary
##==========================================================

cat("\n")
cat("=========================================\n")
cat("Normalization completed successfully\n")
cat("-----------------------------------------\n")
cat("Method : edgeR (TMM)\n")
cat("Transform : limma-voom\n")
cat("-----------------------------------------\n")
cat("Normalized matrix saved to:\n")
cat("results/normalized_expression_matrix.txt\n")
cat("=========================================\n")
