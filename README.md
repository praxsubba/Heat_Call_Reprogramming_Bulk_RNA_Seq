This repository documents a unified, publication-ready bioinformatic workflow to analyze developmental reprogramming in heat-exposed embryos using bulk RNA-sequencing, including differential expression, co-expression network analysis, cell-type deconvolution, and alternative splicing/isoform usage analysis.

***

## Overview

This code base supports four **integrated** RNA-seq analysis streams:

- A **global DESeq2 differential expression pipeline** for bulk RNA-seq.  
- A **targeted hypothalamus DESeq2 pipeline** for an a priori gene set.  
- An **RNA-sequencing co-expression network and deconvolution pipeline**, including WGCNA, module–trait relationships, enrichment, and CIBERSORTx-based deconvolution, which remains unchanged from the original workflow.  
- An **isoform switch and alternative splicing analysis pipeline** based on IsoformSwitchAnalyzeR.

***

## A. Global RNA-seq differential expression (DESeq2)

This component provides a robust, end-to-end DESeq2 analysis for bulk RNA-seq (medial punch) data using Salmon quantifications as input.

### Inputs

**Quantifications**

- `quants/` directory with one subdirectory per sample (e.g., `Sample01_quant/`), each containing `quant.sf` from Salmon.

**Transcript–gene mapping**

- `tx_id.txt`: vector of transcript IDs.  
- `gene_id.txt`: parallel vector of gene IDs used to construct `tx2gene`.

**Sample metadata**

- `sampleInfo.csv` (or `sampleInfo_2.csv`) with at least:
  - `SampleID` (or first column): matches sample identifiers in `quants/`.  
  - `condition`: experimental group (e.g., control vs heat/treatment).  
  - `sex`: biological sex.

### Core script

- `global_rnaseq_deseq2_pipeline.R`  
  This script implements the global DESeq2 workflow and is written to be path-agnostic (no user-specific directories) and safe for repeated, automated runs (timestamped, no-overwrite).

### Analysis steps

**Import and harmonization**

- Builds a transcript-to-gene mapping (`tx2gene`) from `tx_id.txt` and `gene_id.txt`.  
- Aggregates Salmon-level quantifications to gene level with `tximport`.  
- Cleans and matches sample IDs between `quants/` and `sampleInfo`, including handling legacy prefixes and zero-padding where needed.

**DESeq2 modeling**

- Constructs a `DESeqDataSet` with design `~ sex + condition`.  
- Filters low-count genes and fits the treatment model using DESeq2.  
- Performs a likelihood ratio test (LRT) for the condition effect with reduced model `~ sex` for robust assessment of treatment-associated expression changes.

**Log2 fold-change shrinkage and ranking**

- Identifies the appropriate `condition_*` contrast from `resultsNames`.  
- Performs log2 fold-change shrinkage using `apeglm`, generating a ranked table of genes suitable for effect size reporting and downstream interpretation.

**Quality control and visualization**

- Generates a volcano plot of treatment vs control, highlighting a small, well-defined set of labeled genes (top up- and down-regulated).  
- Computes a variance-stabilizing transformation (VST), PCA plot (colored by condition and sex), and a scree plot of principal components.

**Outputs (publication and GEO-ready)**

- Gene-level matrices: raw counts, DESeq2-normalized counts, and VST expression.  
- Sample metadata aligned to matrix columns.  
- Global DESeq2 result tables: unshrunken (full statistics) and shrunken (effect sizes).  
- Optional per-sample CSV files for each data type (raw, normalized, VST), with informative filenames encoding condition, sex, and replicate.

All outputs are written to an `outputs/` directory, with automatic timestamping to prevent overwriting previous runs.

***

## B. Targeted hypothalamus gene-set DESeq2 analysis

This component provides a focused DESeq2 analysis restricted to a **predefined hypothalamus gene list**, allowing rigorous hypothesis-driven testing within an a priori gene set.[3]

### Inputs

In addition to the inputs used for the global pipeline:

**`DeSeq2_hypothalamus_genes.csv`**

- Single-column CSV with header `gene`, containing hypothalamus candidate genes (symbols or IDs).  
- This list defines the targeted analysis universe and is derived independently from the bulk RNA-seq data.[3]

### Core script

- `targeted_hypothalamus_deseq2_pipeline.R`  
  This script reuses the same quantifications and metadata cleaning logic as the global pipeline but restricts analysis to the hypothalamus gene set.[3]

### Analysis steps

**Re-import and metadata harmonization**

- Uses `tximport` to re-import gene-level counts.  
- Cleans and matches sample IDs exactly as in the global pipeline, enforcing consistent handling of `condition` and `sex`.[3]

**Construction of treatment model**

- Builds a `DESeqDataSet` with design `~ sex + condition`.  
- Filters low-count genes prior to restriction to the hypothalamus list.[3]

**Gene-set restriction and testing**

- Restricts the DESeq2 object to genes present both in the dataset and in `DeSeq2_hypothalamus_genes.csv`.  
- Performs DESeq2 with `alpha = 0.05`, yielding unshrunken Wald statistics, followed by `apeglm`-based shrinkage of log2 fold changes.  
- Summarizes the number and direction of significantly altered hypothalamus genes.[3]

**Outputs for figures and supplements**

- Timestamped CSV and RDS files for unshrunken results and shrunken LFCs.  
- Volcano plots for the hypothalamus gene set, using combined FDR and fold-change thresholds.  
- Rank files and combined tables suitable as supplementary resources in the manuscript (unshrunken statistics plus shrunken effect sizes).[3]

This targeted analysis is explicitly labeled as **a priori**, and the hypothalamus gene list is clearly separated from data-driven discovery performed in the global analysis.[3]

***

## C. RNA-sequencing analysis: co-expression network and deconvolution

This section of the pipeline is unchanged from the original workflow and operates on normalized bulk RNA-seq expression matrices generated upstream.

### Inputs

- Raw RNA-seq count data or normalized expression values.  
- Gene length information, if required by specific tools.  
- Sample metadata with experimental traits and covariates.

### Analysis steps

**Setup and data preparation**

- Load count data, gene lengths, and metadata.  
- Harmonize sample IDs and ensure consistent trait encoding.

**Normalization and quality control**

- Apply DESeq2 VST normalization.  
- Perform sample QC and remove outliers as needed.

**Network construction and module detection (WGCNA)**

- Construct a weighted gene co-expression network.  
- Detect gene modules (clusters of co-expressed genes).

**Module–trait relationships**

- Correlate module eigengenes with experimental traits (e.g., condition, sex, phenotypes).  
- Identify modules associated with exposure or developmental outcomes.

**Gene ontology enrichment**

- Perform GO or pathway enrichment on modules of interest.  
- Summarize biological processes and pathways represented by key modules.

**Hub gene identification and validation**

- Identify hub genes within selected modules based on connectivity and module membership.  
- Optionally validate hub genes using external datasets or literature.

**Network visualization**

- Visualize specific modules (e.g., a “green” module) using network plotting tools or Cytoscape-compatible exports.

**Cell-type enrichment analysis**

- Test for enrichment of module genes in cell-type–specific signatures to infer cellular context.

**Single-cell RNA-seq deconvolution signature preparation**

- Generate bulk deconvolution signatures from single-cell reference data.  
- Format signatures for CIBERSORTx or similar deconvolution tools.

**Seurat marker gene identification (optional)**

- Identify marker genes from single-cell or single-nucleus datasets using Seurat, when applicable.

**CIBERSORTx deconvolution (SLURM)**

- Provide a SLURM batch script and input templates for running CIBERSORTx on HPC.  
- Estimate cell-type proportions in bulk RNA-seq samples and integrate these with module-level results.

***

## D. RNA-seq differential expression framework

This section summarizes the general differential expression framework used by the global and targeted RNA-seq analyses.

### Setup and data preparation

**Inputs**

- Salmon quantification (or equivalent) for all samples.  
- Sample metadata with experimental and covariate information.  
- Transcript-to-gene mapping files.  
- Optional tissue-specific gene lists (e.g., hypothalamus genes) used in targeted analyses.[3]

**Tasks**

- Load required R packages.  
- Configure the working directory without user-specific hard-coding.  
- Assemble metadata and mapping files in the expected format.

### Data import and preprocessing

**Tools:** `tximport`, `DESeq2`, `limma`, `qusage`, `dplyr`, `tibble`, `readr`, `ggplot2`, `apeglm`, `ggsci`.

**Tasks**

- Import RNA-seq quantifications and combine with metadata.  
- Construct inputs for DESeq2 and, where appropriate, for downstream model-based analyses.  
- Apply optional tissue-specific filtering (e.g., restricting to hypothalamus gene sets) when performing targeted analyses.[3]

### Differential expression analysis

- Construct DESeq2 objects and filter low-count genes.  
- Fit appropriate models (e.g., `~ sex + condition`).  
- Compute differential expression contrasts and shrink log2 fold changes for robust ranking and interpretation.

### Quality control and visualization

- Visualize sample relationships via PCA and variance-explained plots.  
- Generate volcano plots, rank-based visualizations, and module-/pathway-level figures suitable for inclusion in the main text and supplementary materials.  
- Save all result tables and figures in stable, publication-ready formats (CSV/TSV, PDF/PNG).[3]

***

## E. Isoform switch and alternative splicing analysis

This component implements an isoform-level analysis to detect differential isoform usage and alternative splicing events, with a focus on biologically relevant genes in the developmental heat exposure paradigm.

### Inputs

**Quantifications for isoform-level analysis**

- `quants/` directory used for Salmon quantification, organized with one subdirectory per sample containing Salmon outputs compatible with IsoformSwitchAnalyzeR.

**Sample metadata**

- `sampleInfo.csv` with at least:
  - `SampleID` (or equivalent) matching the samples included in the isoform-level design.  
  - `Comments` or `condition` to encode experimental group.  
  - `sex` for sex as a covariate.

**Reference annotation and sequences**

- Gene annotation (GTF/GFF) with transcript-level features.  
- Transcript FASTA file containing nucleotide sequences for all isoforms.

### Core script

- `splicing_isoform_switch_pipeline.R`  
  This script performs isoform import, design matrix construction, IsoformSwitchAnalyzeR preprocessing, isoform switch testing, splicing characterization, and export of summary tables and gene-level plots.

### Analysis steps

**Import isoform-level expression**

- Uses `importIsoformExpression()` (IsoformSwitchAnalyzeR) to import isoform-level counts and abundances from the Salmon quantification directory.  
- Constructs a `designMatrix` encoding sample IDs, condition, and sex.

**Creation and filtering of the switchAnalyzeRlist**

- Builds a `switchAnalyzeRlist` using isoform counts, abundances, annotation (GTF), and transcript FASTA.  
- Applies `preFilter()` to remove lowly expressed genes and isoforms and to drop single-isoform genes when appropriate for isoform switching tests.

**Isoform switch testing and splicing analysis**

- Runs `isoformSwitchTestDEXSeq()` to identify genes with significant isoform usage changes (differential isoform fraction) between conditions.  
- Extracts summaries of global splicing and alternative splicing categories (e.g., intron retention) using `extractSplicingSummary()` and related functions.  
- Performs sequence extraction and ORF/splicing characterization for detected switches, enabling interpretation at the protein and domain level.

**Gene-focused exploration and outputs**

- Generates switch plots and isoform usage plots for selected genes of interest (e.g., TPM1, TLK2, candidate regulatory genes).  
- Exports a ranked table of isoform switches (e.g., isoforms with `isoform_switch_q_value < 0.01`) and a global summary of splicing events suitable for manuscript tables and supplements.

***

