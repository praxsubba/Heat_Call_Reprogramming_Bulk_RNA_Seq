# Minimal data package for RNA‑seq analysis

This repository contains a minimal set of files required to document and reproduce the core bulk RNA‑seq analyses described in the associated manuscript, limited to the hypothalamus gene list, sample metadata, environment information, and a high‑level description.

---

## Files in this repository

### 1. `README.md`

- High-level description of the repository contents and their role in the RNA‑seq analysis workflow.  
- Intended to be cited in the manuscript as the entry point for code and metadata.

### 2. `DeSeq2_hypothalamus_genes.csv`

- One-column CSV file with header `gene`.  
- Contains the a priori hypothalamus candidate gene list used to define the targeted DESeq2 analysis universe.  
- These genes are derived independently of the bulk RNA‑seq data and are used to restrict the differential expression analysis in the targeted hypothalamus pipeline.

### 3. `sample_metadata.csv`

- Sample-level metadata table for the bulk RNA‑seq dataset.  
- Each row corresponds to a sample, and columns typically include:
  - `sample` or `SampleID`: unique sample identifier matching the expression matrices.  
  - `condition`: experimental group (e.g., control vs treatment).  
  - `sex`: biological sex of the embryo.  
  - Additional covariates used in the DESeq2 models (if applicable).
- This file is suitable for deposition as the “sample metadata” table in GEO/SRA and for direct use with the analysis scripts.

### 4. `sessionInfo_global.txt`

- Text export of `sessionInfo()` from the global RNA‑seq DESeq2 analysis session.  
- Documents:
  - R version and platform.  
  - Versions of all loaded packages (e.g., DESeq2, tximport, ggplot2).  
  - Locale and other session-level settings.[
- This file provides the environment provenance required by many journals for computational reproducibility.

---
