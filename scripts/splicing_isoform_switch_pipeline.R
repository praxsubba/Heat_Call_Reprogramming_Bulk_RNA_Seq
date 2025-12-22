#!/usr/bin/env Rscript

# Isoform switch and alternative splicing analysis pipeline
# Requirements:
#   - R >= 4.x with BiocManager
#   - IsoformSwitchAnalyzeR, matrixStats, Biostrings
#   - Salmon isoform-level quantifications per sample in ./quants/
#   - sampleInfo.csv with SampleID, condition/Comments, sex
#   - Transcript annotation (GTF/GFF) and transcript FASTA

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("IsoformSwitchAnalyzeR", "matrixStats", "Biostrings"),
                       ask = FALSE, update = FALSE)
  library(IsoformSwitchAnalyzeR)
  library(matrixStats)
  library(Biostrings)
})

# --------------------- User-configurable paths ---------------------
# Working directory should contain:
#   - quants/           : Salmon per-sample directories
#   - sampleInfo.csv    : metadata
#   - annotation.gtf.gz : transcript / exon annotation
#   - transcripts.fa.gz : transcript nucleotide FASTA

work_dir        <- "."
quants_dir      <- file.path(work_dir, "quants")
metadata_file   <- file.path(work_dir, "sampleInfo.csv")
gtf_file        <- file.path(work_dir, "annotation.gtf.gz")
tx_fasta_file   <- file.path(work_dir, "transcripts.fa.gz")
results_dir     <- file.path(work_dir, "isoform_switch_results")

if (!dir.exists(quants_dir)) stop("Missing quants/ directory.")
if (!file.exists(metadata_file)) stop("Missing sampleInfo.csv.")
if (!file.exists(gtf_file)) stop("Missing transcript annotation GTF (annotation.gtf.gz).")
if (!file.exists(tx_fasta_file)) stop("Missing transcript FASTA (transcripts.fa.gz).")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------- Import isoform-level expression ---------------------
message("Importing isoform-level expression from Salmon output...")

salmon_quant <- importIsoformExpression(
  parentDir            = quants_dir,
  addIsofomIdAsColumn  = TRUE
)

# --------------------- Design matrix ---------------------
meta <- read.csv(metadata_file, header = TRUE, stringsAsFactors = FALSE)

# Expect a column that maps to Salmon directory names; adjust this logic as needed
if (!"SampleID" %in% colnames(meta)) {
  stop("sampleInfo.csv must contain a 'SampleID' column.")
}

# Example: Salmon directories named "<SampleID>_quant"
sample_ids <- paste0(meta$SampleID, "_quant")

design <- data.frame(
  sampleID  = sample_ids,
  condition = ifelse("condition" %in% colnames(meta), meta$condition,
                     ifelse("Comments" %in% colnames(meta), meta$Comments, NA)),
  sex       = if ("sex" %in% colnames(meta)) meta$sex else NA,
  stringsAsFactors = FALSE
)

if (any(is.na(design$condition))) {
  stop("Condition column not found or contains NA; ensure 'condition' or 'Comments' exists.")
}

# --------------------- Create switchAnalyzeRlist ---------------------
message("Creating switchAnalyzeRlist...")

switch_list <- importRdata(
  isoformCountMatrix   = salmon_quant$counts,
  isoformRepExpression = salmon_quant$abundance,
  designMatrix         = design,
  isoformExonAnnoation = gtf_file,
  isoformNtFasta       = tx_fasta_file,
  showProgress         = FALSE
)

# --------------------- Optional: sanity check for genes of interest ---------------------
# Example: check presence and isoform structure for a specific gene (e.g., NR3C1)
check_gene <- function(gene_name, fasta_file, switch_list_object) {
  tx <- readDNAStringSet(fasta_file)
  gene_names <- sapply(names(tx), function(x) {
    parts <- strsplit(x, " ")[[1]]
    gene_part <- parts[grep("gene=", parts)]
    if (length(gene_part) > 0) {
      gsub("gene=", "", gene_part)
    } else {
      NA
    }
  })
  gene_names <- unique(gene_names[!is.na(gene_names)])

  if (!(gene_name %in% gene_names)) {
    message(gene_name, " not found in transcript FASTA.")
    similar <- gene_names[grep(gene_name, gene_names, ignore.case = TRUE)]
    if (length(similar)) {
      message("Similar gene names: ", paste(similar, collapse = ", "))
    }
    return(invisible(NULL))
  }

  # Extract isoform IDs for gene_name
  tx_gene <- tx[grep(paste0("gene=", gene_name), names(tx))]
  iso_ids <- names(tx_gene)
  message("Number of transcripts for ", gene_name, ": ", length(iso_ids))

  # Map to expression matrix
  mat <- as.matrix(switch_list_object$isoformCountMatrix[, -1])
  rownames(mat) <- switch_list_object$isoformCountMatrix$isoform_id
  matching_ids <- intersect(iso_ids, rownames(mat))

  if (!length(matching_ids)) {
    message("No matching ", gene_name, " isoforms found in count matrix.")
    return(invisible(NULL))
  }

  expr <- mat[matching_ids, , drop = FALSE]
  mean_expr <- rowMeans(expr)
  gene_expr <- sum(mean_expr)
  message(gene_name, " isoform mean expression:"); print(mean_expr)
  message("Overall ", gene_name, " expression: ", gene_expr)

  invisible(list(mean_isoform_expr = mean_expr, total_expr = gene_expr))
}

# Example call (comment out or change gene_name as needed)
# check_gene("NR3C1", tx_fasta_file, switch_list)

# --------------------- Filtering ---------------------
message("Filtering lowly expressed genes and isoforms...")

switch_list_filt <- preFilter(
  switchAnalyzeRlist       = switch_list,
  geneExpressionCutoff     = 1,
  isoformExpressionCutoff  = 0,
  removeSingleIsoformGenes = TRUE
)

# --------------------- Isoform switch testing ---------------------
message("Testing for isoform switches with DEXSeq...")

switch_list_analyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist     = switch_list_filt,
  reduceToSwitchingGenes = TRUE,
  alpha                  = 0.05,
  dIFcutoff              = 0.05
)

# Summary of switching features
switch_summary <- extractSwitchSummary(switch_list_analyzed)
print(switch_summary)

saveRDS(switch_list_analyzed,
        file = file.path(results_dir, "isoformSwitchTestDEXSeq_results.rds"))

# --------------------- Sequence and splicing analysis ---------------------
message("Extracting sequences and performing splicing analysis...")

switch_list_analyzed <- extractSequence(
  switch_list_analyzed,
  pathToOutput = results_dir,
  writeToFile  = FALSE
)

switch_list_analyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = switch_list_analyzed,
  quiet              = TRUE
)

# Global splicing summaries
splicing_summary <- extractSplicingSummary(
  switch_list_analyzed,
  asFractionTotal = FALSE,
  plotGenes       = FALSE
)
print(splicing_summary)

splicing_enrichment <- extractSplicingEnrichment(
  switch_list_analyzed,
  splicingToAnalyze = "all",
  returnResult      = TRUE,
  returnSummary     = TRUE
)
print(splicing_enrichment)

# --------------------- Gene-level plots (examples) ---------------------
# Replace gene symbols below with genes of interest in your study.
genes_of_interest <- c("TPM1", "TLK2")

pdf(file.path(results_dir, "switch_plots_example_genes.pdf"), width = 7, height = 5)
for (g in genes_of_interest) {
  try(switchPlot(switch_list_analyzed, gene = g), silent = TRUE)
}
dev.off()

pdf(file.path(results_dir, "isoform_usage_plots_example_genes.pdf"), width = 7, height = 5)
for (g in genes_of_interest) {
  try(switchPlotIsoUsage(switch_list_analyzed, gene = g), silent = TRUE)
}
dev.off()

# --------------------- Export ranked isoform switches ---------------------
isoform_features <- switch_list_analyzed$isoformFeatures
isoform_features_sorted <- isoform_features[order(isoform_features$isoform_switch_q_value), ]

sig_isoforms <- subset(
  isoform_features_sorted,
  !is.na(isoform_switch_q_value) & isoform_switch_q_value < 0.01
)

sig_isoforms_export <- sig_isoforms[, c("isoform_id", "gene_id", "isoform_switch_q_value")]
write.csv(
  sig_isoforms_export,
  file = file.path(results_dir, "isoform_usage_results_qval_lt_0_01.csv"),
  row.names = FALSE
)

# Save splicing summaries for manuscript
write.csv(
  splicing_summary,
  file = file.path(results_dir, "splicing_summary.csv"),
  row.names = FALSE
)
write.csv(
  splicing_enrichment,
  file = file.path(results_dir, "splicing_enrichment.csv"),
  row.names = FALSE
)

# Save complete R environment if desired
save.image(file = file.path(results_dir, "Isoform_Switch_Analysis.RData"))

sessionInfo()
