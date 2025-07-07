# RNA-seq Differential Expression Analysis with DESeq2

## Overview

This project provides a comprehensive R-based pipeline for **RNA-seq differential gene expression analysis** using **DESeq2**, specifically developed for liver tissue transcriptomic data from an animal trial investigating **vitamin D and phosphorus metabolism**. The workflow includes raw count import, annotation, metadata handling, normalization, filtering, differential testing, and outlier diagnostics.

---

## Dataset and Tools

- **Input Files**:
  - `Liver_AFBI-UVB_pheno_subset.csv` – Phenotypic metadata file
  - `annotated.ensembl.gene.name.txt` – Gene annotation file
  - `*.GCount.txt` – Raw HTSeq-count output per sample

- **Tools & Libraries**:
  - R (v4.0 or later)
  - Required packages:
    - `DESeq2`, `qvalue`, `affy`, `genefilter`, `arrayQualityMetrics`, `GMD`

---

## Goals

- Import and clean **HTSeq-count RNA-seq data**
- Link with **phenotype metadata**
- Perform **differential expression analysis** using DESeq2
- Apply **Cook's distance filtering** and **multiple testing correction**
- Annotate results with **gene names**
- Export a publication-ready result table of **significantly differentially expressed genes**

---

## Repository Structure

```
RNAseq_DESeq2_Analysis/
├── annotated.ensembl.gene.name.txt         # Gene annotation mapping file
├── Liver_AFBI-UVB_pheno_subset.csv         # Sample metadata
├── htseq.gene.count/                       # Folder with HTSeq count files (*.GCount.txt)
├── R/
│   ├── RNAseq_DESeq2_Script.R              # Main R analysis script
│   ├── Liver_Group_None_anno_DESeq2.csv    # Final DE result output (example)
```

---

## Reproducibility Instructions

1. **Set the working directory**:
```r
setwd("Z:/FBNProjectGroups/I3-NGS-DexaPhos/Maruf_VitaminD/RNA-seq/DESeq2")
```

2. **Load required R libraries**:
```r
library(DESeq2)
library(qvalue)
library(affy)
library(genefilter)
library(arrayQualityMetrics)
library(GMD)
```

3. **Customize script parameters**:
```r
tissue <- c("Liv")
effect <- c("Group")
outlier <- c("None")
```

4. **Run the script**:
```r
source("RNAseq_DESeq2_Script.R")
```

---

## Key Features

- **Automatic raw count import**:
  - Imports all `*.GCount.txt` files from a given directory
  - Combines them into a single count matrix

- **Phenotype and gene mapping integration**:
  - Merges phenotype data by sample ID
  - Merges gene annotation by Ensembl ID

- **Differential expression analysis**:
  - Uses Wald test with customizable model formula (`~ Group`)
  - Filters low-expression genes (≥50 counts in ≥4 samples)
  - Applies Cook’s distance filtering for outlier control
  - Adjusts p-values using Benjamini-Hochberg (FDR)

- **Final output**:
  - Annotated DE gene list sorted by p-value and adjusted p-value
  - Output CSV: e.g., `Liver_Group_None_anno_DESeq2.csv`

---

## Output Example

- `Liver_Group_None_anno_DESeq2.csv` includes:
  - `SSC_ID`: Ensembl Gene ID
  - `log2FoldChange`: Effect size
  - `pvalue`: Raw p-value
  - `padj`: Adjusted p-value (FDR)
  - `gene`: Gene symbol

---

## Diagnostics

- **Dispersion plot**: `plotDispEsts(dds)`
- **Cook’s distance filtering**:
  - Max Cook’s distance per gene
  - Outlier identification and filtering based on 99.99th percentile
- **Sanity checks**:
  - Match of DESeq2 rownames across `dds` and `res`

---

## Purpose

This RNA-seq pipeline is tailored for hypothesis testing of **dietary or genetic effects** on gene expression in **livestock models**, particularly under conditions affecting **mineral metabolism** and **vitamin D regulation**.

---

## Author

**Maruf Hasan**  
Interests: Transcriptomics | Vitamin D | Mineral metabolism | Livestock models 
