############################################################
# Project: TCGA-RNAseq-Analysis
# Script : 04_differential_expression.R
# Author : Zeinab Rohani
#
# Description:
# Differential gene expression analysis using
# limma and empirical Bayes statistics.
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
## Load Normalized Expression Matrix
##==========================================================

expression_matrix <- read.delim(
  "results/normalized_expression_matrix.txt",
  row.names = 1,
  check.names = FALSE
)

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
## Create Design Matrix
##==========================================================

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

##==========================================================
## Fit Linear Model
##==========================================================

fit <- lmFit(
  expression_matrix,
  design
)

##==========================================================
## Define Comparison
##==========================================================

contrast_matrix <- makeContrasts(
  Cancer - Normal,
  levels = design
)

##==========================================================
## Fit Contrast
##==========================================================

fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)

##==========================================================
## Empirical Bayes Moderation
##==========================================================

fit2 <- eBayes(fit2)

##==========================================================
## Extract Differentially Expressed Genes
##==========================================================

deg_results <- topTable(
  fit2,
  number = Inf,
  sort.by = "P"
)

##==========================================================
## Add Gene IDs
##==========================================================

deg_results$Ensembl_ID <- rownames(deg_results)

deg_results <- deg_results[
  ,
  c(
    "Ensembl_ID",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
]

##==========================================================
## Export Results
##==========================================================

write.table(
  deg_results,
  file = "results/differential_expression_results.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

##==========================================================
## Significant DEGs
##==========================================================

significant_deg <- subset(
  deg_results,
  abs(logFC) >= 1 &
  adj.P.Val < 0.05
)

write.table(
  significant_deg,
  file = "results/significant_DEGs.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

##==========================================================
## Summary
##==========================================================

cat("\n")
cat("=========================================\n")
cat("Differential Expression Analysis Complete\n")
cat("-----------------------------------------\n")
cat("Total genes tested :", nrow(deg_results), "\n")
cat("Significant DEGs   :", nrow(significant_deg), "\n")
cat("-----------------------------------------\n")
cat("Results exported to results/\n")
cat("=========================================\n")
