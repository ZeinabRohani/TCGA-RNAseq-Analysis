############################################################
# Project: TCGA-RNAseq-Analysis
# Script : 01_download_TCGA.R
# Author : Zeinab Rohani
# Description:
# Download TCGA-COAD RNA-seq count data and prepare a
# SummarizedExperiment object using TCGAbiolinks.
############################################################

##==========================================================
## Load Required Packages
##==========================================================

library(TCGAbiolinks)
library(SummarizedExperiment)

##==========================================================
## Set Working Directory
##==========================================================

setwd("D:/TCGA/")

##==========================================================
## Create Project Folders (if they do not exist)
##==========================================================

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

##==========================================================
## Query TCGA Database
##==========================================================

query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq",
  legacy = FALSE
)

##==========================================================
## Download Data
##==========================================================

GDCdownload(
  query,
  method = "api",
  directory = "data/",
  files.per.chunk = 10
)

##==========================================================
## Prepare SummarizedExperiment Object
##==========================================================

tcga_data <- GDCprepare(
  query,
  directory = "data/",
  summarizedExperiment = TRUE
)

##==========================================================
## Save Object
##==========================================================

saveRDS(
  tcga_data,
  file = "data/tcga_coad_summarized_experiment.rds"
)

##==========================================================
## Summary
##==========================================================

cat("\n")
cat("=========================================\n")
cat("TCGA COAD RNA-seq data downloaded\n")
cat("SummarizedExperiment object created\n")
cat("Saved to:\n")
cat("data/tcga_coad_summarized_experiment.rds\n")
cat("=========================================\n")
