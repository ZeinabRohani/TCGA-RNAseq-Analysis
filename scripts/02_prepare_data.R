############################################################
# Project: TCGA-RNAseq-Analysis
# Script : 02_prepare_data.R
# Author : Zeinab Rohani
# Description:
# Extract count matrix, clinical metadata,
# and gene annotation from the
# SummarizedExperiment object.
############################################################

##==========================================================
## Load Required Packages
##==========================================================

library(SummarizedExperiment)

##==========================================================
## Set Working Directory
##==========================================================

setwd("D:/TCGA/")

##==========================================================
## Create Results Folder
##==========================================================

dir.create("results", showWarnings = FALSE)

##==========================================================
## Load SummarizedExperiment Object
##==========================================================

tcga_data <- readRDS("data/tcga_coad_summarized_experiment.rds")

##==========================================================
## Extract Raw Count Matrix
##==========================================================

count_matrix <- assay(tcga_data)

write.table(
  count_matrix,
  file = "results/raw_counts.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

##==========================================================
## Extract Clinical Metadata
##==========================================================

clinical_data <- as.data.frame(colData(tcga_data))

write.table(
  clinical_data,
  file = "results/clinical_metadata.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

##==========================================================
## Extract Gene Annotation
##==========================================================

gene_annotation <- as.data.frame(rowRanges(tcga_data))

write.table(
  gene_annotation,
  file = "results/gene_annotation.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

##==========================================================
## Print Summary
##==========================================================

cat("\n")
cat("=========================================\n")
cat("Data preparation completed successfully\n")
cat("-----------------------------------------\n")
cat("Raw counts exported\n")
cat("Clinical metadata exported\n")
cat("Gene annotation exported\n")
cat("-----------------------------------------\n")
cat("Files saved in: results/\n")
cat("=========================================\n")
