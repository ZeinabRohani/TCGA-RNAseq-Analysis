# TCGA-RNAseq-Analysis

> **A reproducible RNA-seq differential gene expression analysis pipeline for TCGA colorectal cancer (TCGA-COAD) using R and Bioconductor.**

![R](https://img.shields.io/badge/R-4.3+-276DC3?logo=r)
![Bioconductor](https://img.shields.io/badge/Bioconductor-3.x-green)
![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux-blue)
![License](https://img.shields.io/badge/License-MIT-yellow)

---

# Overview

This repository provides a reproducible workflow for analyzing RNA sequencing (RNA-seq) data from **The Cancer Genome Atlas (TCGA)** using **R** and **Bioconductor**.

The pipeline demonstrates the complete process of downloading transcriptomic data, preprocessing raw counts, filtering lowly expressed genes, normalization, differential gene expression analysis, and gene annotation.

The workflow is implemented using widely adopted Bioconductor packages including **TCGAbiolinks**, **edgeR**, **limma**, and **SummarizedExperiment**, following best practices for bulk RNA-seq analysis.

This project was developed as part of my computational biology portfolio and reflects my experience in cancer transcriptomics and bioinformatics.

---

# Dataset

- **Database:** The Cancer Genome Atlas (TCGA)
- **Cancer Type:** Colon Adenocarcinoma (TCGA-COAD)
- **Data Type:** Gene Expression Quantification
- **Workflow:** STAR - Counts / HTSeq - Counts
- **Platform:** RNA Sequencing

---

# Objectives

The objectives of this project are to:

- Download RNA-seq count data from TCGA
- Retrieve clinical metadata
- Extract gene annotation
- Filter lowly expressed genes
- Normalize count data
- Perform differential gene expression analysis
- Annotate Ensembl gene IDs
- Export reproducible analysis results

---

# Workflow

```text
                    TCGA COAD RNA-seq
                           │
                           ▼
                 Download Data (TCGAbiolinks)
                           │
                           ▼
            Prepare SummarizedExperiment Object
                           │
            ┌──────────────┴──────────────┐
            ▼                             ▼
      Clinical Data                Gene Annotation
            │                             │
            └──────────────┬──────────────┘
                           ▼
                Raw Count Matrix Extraction
                           │
                           ▼
               Filter Low-Expression Genes
                           │
                           ▼
              TMM Normalization (edgeR)
                           │
                           ▼
               Voom Transformation (limma)
                           │
                           ▼
                Linear Model Construction
                           │
                           ▼
              Empirical Bayes Moderation
                           │
                           ▼
          Differential Gene Expression Analysis
                           │
                           ▼
               Gene Symbol Annotation
                           │
                           ▼
                    Export Results
```

---

# Repository Structure

```text
TCGA-RNAseq-Analysis/

├── data/
│
├── results/
│   ├── expression matrix
│   ├── normalized matrix
│   ├── differential expression
│   ├── clinical data
│   └── gene annotation
│
├── scripts/
│   ├── Download_TCGA.R
│   ├── Normalization.R
│   ├── Differential_Expression.R
│   └── Annotation.R
│
├── README.md
└── LICENSE
```

---

# Software

- R (≥4.3)
- RStudio

---

# Main Packages

| Package | Purpose |
|----------|---------|
| TCGAbiolinks | Download TCGA RNA-seq data |
| edgeR | Filtering and TMM normalization |
| limma | Differential expression analysis |
| SummarizedExperiment | Expression data management |

---

# Analysis Pipeline

## 1. Data Acquisition

RNA-seq count data were downloaded directly from TCGA using **TCGAbiolinks**.

---

## 2. Data Preparation

The workflow extracts:

- Raw count matrix
- Clinical metadata
- Gene annotation

from the SummarizedExperiment object.

---

## 3. Gene Filtering

Lowly expressed genes are removed using

- `filterByExpr()`

to improve statistical power.

---

## 4. Normalization

Library size differences are corrected using

- TMM normalization (`calcNormFactors()`)

followed by

- `voom()` transformation

for downstream linear modeling.

---

## 5. Differential Expression Analysis

Differentially expressed genes are identified using

- Linear models (`lmFit`)
- Contrast matrix
- Empirical Bayes moderation (`eBayes`)
- `topTable()`

---

## 6. Gene Annotation

Ensembl gene identifiers are mapped to official HGNC gene symbols.

---

# Outputs

The pipeline generates:

- Raw count matrix
- Clinical metadata
- Gene annotation table
- Normalized expression matrix
- Differentially expressed gene table

---

# Skills Demonstrated

- Cancer Transcriptomics
- RNA-seq Analysis
- TCGA Data Retrieval
- Differential Gene Expression
- Statistical Modeling
- Empirical Bayes Statistics
- Gene Annotation
- edgeR
- limma
- TCGAbiolinks
- Bioconductor
- R Programming
- Reproducible Bioinformatics

---

# Related Publications

The computational approaches demonstrated in this repository are closely related to my research in colorectal cancer transcriptomics.

### Unlocking the potential of *Escherichia coli* K-12: A novel approach for malignancy reduction in colorectal cancer through gene expression modulation

https://scholar.google.com/citations?view_op=view_citation&hl=en&user=cRFVNUQAAAAJ&citation_for_view=cRFVNUQAAAAJ:u5HHmVD_uO8C

---

### In Silico Identification of Phenylacetaldehyde's Impact on Downregulating Malignancy-Associated Genes: A Promising Therapeutic Approach Against Colorectal Cancer

https://scholar.google.com/citations?view_op=view_citation&hl=en&user=cRFVNUQAAAAJ&citation_for_view=cRFVNUQAAAAJ:d1gkVwhDpl0C

---

# Citation

If you use this repository in your research, please cite:

- The Cancer Genome Atlas (TCGA)
- TCGAbiolinks
- edgeR
- limma

---

# Author

**Zeinab Rohani**

**Research Interests**

- Bioinformatics
- Computational Biology
- Cancer Genomics
- Transcriptomics
- Functional Genomics
- RNA-seq Analysis
