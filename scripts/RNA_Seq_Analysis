# Global RNA-seq differential expression analysis
# Requirements:
# - Working directory contains:
#     quants/            # per-sample Salmon quant directories
#     tx_id.txt          # transcript IDs from reference
#     gene_id.txt        # corresponding gene IDs
#     sampleInfo.csv     # or sampleInfo_2.csv, with SampleID, condition, sex
# - Outputs written to ./outputs/ with timestamped, no-overwrite behavior.

# -------------------- Setup --------------------
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("tximport", "DESeq2", "apeglm"), ask = FALSE, update = FALSE)
  install.packages(c("dplyr", "ggplot2", "tidyr", "tibble",
                     "readr", "forcats", "ggsci", "jsonlite", "ggrepel"),
                   quiet = TRUE)
  library(dplyr)
  library(tximport)
  library(DESeq2)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(apeglm)
  library(ggsci)
  library(tidyr)
  library(forcats)
  library(ggrepel)
})

# -------------------- Helper functions --------------------
dir_create_if_missing <- function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

timestamped <- function(path) {
  if (!file.exists(path)) return(path)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  file.path(
    dirname(path),
    paste0(tools::file_path_sans_ext(basename(path)), "_", ts, ".", tools::file_ext(path))
  )
}

save_csv_no_overwrite <- function(df, path) {
  tgt <- if (!file.exists(path)) path else timestamped(path)
  write.csv(df, tgt, row.names = FALSE)
  message("Saved: ", tgt)
  invisible(tgt)
}

save_rds_no_overwrite <- function(obj, path) {
  tgt <- if (!file.exists(path)) path else timestamped(path)
  saveRDS(obj, tgt)
  message("Saved: ", tgt)
  invisible(tgt)
}

ggsave_no_overwrite <- function(filename, plot, width = 7.5, height = 5.5, dpi = 300) {
  tgt <- if (!file.exists(filename)) filename else {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    paste0(tools::file_path_sans_ext(filename), "_", ts, ".", tools::file_ext(filename))
  }
  ggsave(tgt, plot, width = width, height = height, dpi = dpi)
  message("Saved: ", tgt)
  invisible(tgt)
}

# -------------------- Paths and validation --------------------
outputs_dir <- "outputs"
dir_create_if_missing(outputs_dir)

if (!all(file.exists(c("tx_id.txt", "gene_id.txt")))) {
  stop("Missing tx_id.txt or gene_id.txt in working directory.")
}
if (!dir.exists("quants")) {
  stop("Missing quants/ directory next to this script.")
}

# -------------------- tx2gene mapping --------------------
tx_raw   <- read.csv("tx_id.txt", header = FALSE)$V1
gene_raw <- read.csv("gene_id.txt", header = FALSE)$V1
tx2gene  <- data.frame(
  TXNAME = gsub("\\.\\d+$", "", tx_raw),
  GENEID = gene_raw,
  stringsAsFactors = FALSE
)
message("tx2gene ready: n=", nrow(tx2gene))

# -------------------- Import Salmon quantifications --------------------
quants_dir  <- "quants"
quant_dirs  <- list.files(quants_dir)
if (length(quant_dirs) == 0) stop("No sample subdirectories found in quants/")

files <- file.path(quants_dir, quant_dirs, "quant.sf")
if (!all(file.exists(files))) {
  missing <- files[!file.exists(files)]
  stop("Missing quant.sf files, examples: ", paste(head(missing, 3), collapse = ", "))
}

sampleID <- gsub("_quant$", "", quant_dirs)
names(files) <- sampleID

quant_head <- readr::read_tsv(files[1], n_max = 6, col_names = TRUE, show_col_types = FALSE)
message("Example quant IDs: ", paste(head(quant_head$Name, 5), collapse = ", "))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
message("tximport complete: counts dims = ", paste(dim(txi$counts), collapse = " x "))

# -------------------- Metadata --------------------
csv_file <- if (file.exists("sampleInfo.csv")) {
  "sampleInfo.csv"
} else if (file.exists("sampleInfo_2.csv")) {
  "sampleInfo_2.csv"
} else {
  stop("Missing sampleInfo.csv or sampleInfo_2.csv")
}

treat <- read.csv(csv_file, header = TRUE)

id_col_name <- ifelse("SampleID" %in% names(treat), "SampleID", names(treat)[1])
condition_col_name <- ifelse(
  "condition" %in% names(treat), "condition",
  ifelse("Comments" %in% names(treat), "Comments",
         stop("No 'condition' or 'Comments' column in metadata"))
)
sex_col_name <- "sex"

id_col <- treat[[id_col_name]]

treat_ids_clean <- as.character(id_col)
treat_ids_clean <- gsub("^[S|s]?", "", treat_ids_clean)
treat_ids_clean <- gsub("[_ ]+$", "", treat_ids_clean)
treat_ids_clean <- gsub("^PREFIX-", "", treat_ids_clean)  # study-specific prefix (anonymized)
treat_ids_clean <- trimws(treat_ids_clean)
treat_ids_clean <- ifelse(
  nchar(gsub("[^0-9]", "", treat_ids_clean)) == 1 &
    grepl("^[0-9]+$", treat_ids_clean),
  paste0("0", treat_ids_clean),
  treat_ids_clean
)

samples <- data.frame(
  sampleID  = sampleID,
  condition = treat[[condition_col_name]][match(sampleID, treat_ids_clean)],
  sex       = treat[[sex_col_name]][match(sampleID, treat_ids_clean)],
  stringsAsFactors = FALSE
)
samples$sex_condition <- factor(paste0(samples$sex, samples$condition))

na_count <- sum(is.na(samples$condition) | is.na(samples$sex))
if (na_count > 0) {
  mismatches <- sampleID[is.na(samples$condition) | is.na(samples$sex)]
  stop("Metadata mismatches (", na_count, ") e.g., ",
       paste(head(mismatches, 5), collapse = ", "),
       ". Adjust cleaning or fix the CSV.")
}
message("Metadata matched: n=", nrow(samples))

# -------------------- DESeq2: treatment model --------------------
dds_treat_path <- file.path(outputs_dir, "dds_treatment.rds")
if (file.exists(dds_treat_path)) {
  dds_treat <- readRDS(dds_treat_path)
  message("Loaded existing dds_treat: ", dds_treat_path)
} else {
  dds_treat <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ sex + condition)
  keep <- rowSums(counts(dds_treat)) >= 10
  dds_treat <- dds_treat[keep, ]
  set.seed(1)
  dds_treat <- DESeq(dds_treat)
  save_rds_no_overwrite(dds_treat, dds_treat_path)
}

# Ensure factors and refit
colData(dds_treat)$sex <- factor(colData(dds_treat)$sex)
if (is.numeric(colData(dds_treat)$condition) || is.integer(colData(dds_treat)$condition)) {
  colData(dds_treat)$condition <- factor(
    colData(dds_treat)$condition,
    levels = c(0, 1),
    labels = c("control", "treatment")
  )
} else {
  cd <- tolower(as.character(colData(dds_treat)$condition))
  cd[cd %in% c("control", "ctrl", "c", "0")] <- "control"
  cd[cd %in% c("treatment", "treated", "t", "1", "heat")] <- "treatment"
  colData(dds_treat)$condition <- factor(cd, levels = c("control", "treatment"))
}
set.seed(1)
dds_treat <- DESeq(dds_treat, quiet = TRUE)

# Likelihood ratio test (treatment effect)
dds_lrt_treat <- DESeq(dds_treat, test = "LRT", reduced = ~ sex, quiet = TRUE)
res_treat <- results(dds_lrt_treat, alpha = 0.05)
res_treat_df <- as.data.frame(res_treat) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)
save_csv_no_overwrite(res_treat_df, file.path(outputs_dir, "DESeq2_results_all_genes.csv"))

# -------------------- LFC shrinkage --------------------
rn <- resultsNames(dds_treat)
message("resultsNames: ", paste(rn, collapse = " | "))

coef_name <- NULL
cand <- rn[grepl("^condition_", rn)]
if (length(cand)) {
  if (any(grepl("treatment_vs_control", cand))) {
    coef_name <- cand[grepl("treatment_vs_control", cand)][1]
  } else if (any(grepl("_treatment_vs_", cand))) {
    coef_name <- cand[grepl("_treatment_vs_", cand)][1]
  } else {
    coef_name <- cand[1]
  }
} else {
  stop("No condition_* coefficients found in resultsNames: ", paste(rn, collapse = ", "))
}
message("Using coef for shrinkage: ", coef_name)

res_global <- results(dds_treat, alpha = 0.05)
res_lfc <- lfcShrink(dds_treat, res = res_global, coef = coef_name, type = "apeglm")
res_lfc_df <- as.data.frame(res_lfc) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)
save_csv_no_overwrite(res_lfc_df, file.path(outputs_dir, "DESeq2_results_shrunken_LFC.csv"))
save_rds_no_overwrite(res_lfc, file.path(outputs_dir, "DESeq2_results_shrunken_LFC.rds"))

# -------------------- Volcano plot --------------------
volcano_df <- res_treat_df %>%
  mutate(
    neg_log_pvalue = -log10(pvalue),
    direction = dplyr::case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Up Sig",
      padj < 0.05 & log2FoldChange < 0 ~ "Down Sig",
      TRUE ~ "Not Sig"
    )
  )

top5_down <- volcano_df %>%
  filter(padj < 0.05, log2FoldChange < 0) %>%
  arrange(padj) %>%
  slice_head(n = 5) %>%
  pull(gene)

top1_up <- volcano_df %>%
  filter(padj < 0.05, log2FoldChange > 0) %>%
  arrange(padj) %>%
  slice_head(n = 1) %>%
  pull(gene)

label_genes <- unique(c(top5_down, top1_up))

volcano_df <- volcano_df %>%
  mutate(label = ifelse(gene %in% label_genes, gene, NA_character_))

p_volcano <- ggplot(volcano_df,
                    aes(x = log2FoldChange, y = neg_log_pvalue, color = direction)) +
  geom_point(alpha = ifelse(volcano_df$direction == "Not Sig", 0.4, 0.8), size = 1.8) +
  scale_color_manual(
    values = c("Down Sig" = "#1f77b4",
               "Up Sig"   = "#d62728",
               "Not Sig"  = "grey60"),
    guide = guide_legend(override.aes = list(alpha = 1, size = 2))
  ) +
  geom_text_repel(
    data = dplyr::filter(volcano_df, !is.na(label)),
    aes(label = label),
    max.overlaps = 10,
    size = 4.0,
    box.padding = 1.0,
    point.padding = 1.0,
    segment.color = "black",
    segment.size = 0.35,
    segment.alpha = 0.6,
    fontface = "bold",
    seed = 42
  ) +
  labs(
    title = "Volcano Plot: Treatment Effect on Gene Expression",
    x = expression(log[2] ~ "Fold Change (treatment vs control)"),
    y = expression(-log[10] ~ "(p-value)")
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

ggsave(
  filename = file.path(outputs_dir, "volcano_treatment_effect.pdf"),
  plot     = p_volcano,
  width    = 7.5,
  height   = 5.5,
  dpi      = 300,
  device   = cairo_pdf
)
ggsave_no_overwrite(
  file.path(outputs_dir, "volcano_treatment_effect.png"),
  p_volcano,
  width = 8.5,
  height = 7,
  dpi = 600
)

# -------------------- PCA and scree --------------------
vsd <- vst(dds_treat, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = c("condition", "sex"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA of RNA-seq Samples") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()
ggsave_no_overwrite(file.path(outputs_dir, "pca_scatter.pdf"), p_pca, width = 7, height = 5, dpi = 300)
ggsave_no_overwrite(file.path(outputs_dir, "pca_scatter.png"), p_pca, width = 7, height = 5, dpi = 300)

pca_obj <- prcomp(t(assay(vsd)))
scree_df <- data.frame(
  PC = 1:length(pca_obj$sdev),
  var_explained = (pca_obj$sdev^2 / sum(pca_obj$sdev^2)) * 100
)
p_scree <- ggplot(scree_df, aes(x = PC, y = var_explained)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_line(color = "red") +
  geom_point(color = "red") +
  labs(
    x = "Principal Component",
    y = "Percent Variance Explained",
    title = "Scree Plot"
  ) +
  theme_minimal()
ggsave_no_overwrite(file.path(outputs_dir, "pca_scree.pdf"), p_scree, width = 7, height = 5, dpi = 300)
ggsave_no_overwrite(file.path(outputs_dir, "pca_scree.png"), p_scree, width = 7, height = 5, dpi = 300)

# -------------------- GEO-ready exports --------------------
# 1) Raw counts
raw_counts <- counts(dds_treat, normalized = FALSE)
raw_counts_df <- as.data.frame(raw_counts) %>%
  tibble::rownames_to_column("gene")
save_csv_no_overwrite(raw_counts_df, file.path(outputs_dir, "gene_counts_raw.csv"))

# 2) Normalized counts
norm_counts <- counts(dds_treat, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts) %>%
  tibble::rownames_to_column("gene")
save_csv_no_overwrite(norm_counts_df, file.path(outputs_dir, "gene_counts_DESeq2_normalized.csv"))

# 3) VST expression
vst_mat <- assay(vsd)
vst_df <- as.data.frame(vst_mat) %>%
  tibble::rownames_to_column("gene")
save_csv_no_overwrite(vst_df, file.path(outputs_dir, "gene_expression_VST.csv"))

# 4) Sample metadata
sample_metadata_df <- as.data.frame(colData(dds_treat))
sample_metadata_df <- tibble::rownames_to_column(sample_metadata_df, var = "sample")
save_csv_no_overwrite(sample_metadata_df, file.path(outputs_dir, "sample_metadata.csv"))

# 5) Per-sample processed files (raw, normalized, VST)
sample_meta <- as.data.frame(colData(dds_treat))
if (!"sampleID" %in% colnames(sample_meta)) {
  sample_meta <- sample_meta %>% tibble::rownames_to_column("sampleID")
}
sample_meta <- sample_meta %>%
  select(sampleID, condition, sex) %>%
  mutate(
    condition_clean = case_when(
      tolower(as.character(condition)) %in% c("control", "ctrl", "c", "0") ~ "Control",
      tolower(as.character(condition)) %in% c("treatment", "treated", "t", "1", "heat") ~ "Treatment",
      TRUE ~ as.character(condition)
    ),
    sex_clean = case_when(
      toupper(as.character(sex)) == "F" ~ "Female",
      toupper(as.character(sex)) == "M" ~ "Male",
      TRUE ~ as.character(sex)
    ),
    sampleID_padded = sprintf("%02d", as.numeric(gsub("\\D", "", sampleID)))
  )

export_per_sample <- function(matrix_df, sample_meta, data_type_label, outputs_dir) {
  gene_col <- matrix_df[, 1, drop = TRUE]
  data_cols <- matrix_df[, -1, drop = FALSE]
  for (col_name in colnames(data_cols)) {
    meta_row <- sample_meta[sample_meta$sampleID == col_name, ]
    if (nrow(meta_row) == 0) {
      warning("No metadata found for sample: ", col_name, "; skipping.")
      next
    }
    filename <- paste0(
      meta_row$condition_clean, "_",
      meta_row$sex_clean, "_",
      "Rep_", meta_row$sampleID_padded, "_",
      data_type_label, ".csv"
    )
    sample_df <- data.frame(
      gene  = gene_col,
      value = data_cols[[col_name]],
      stringsAsFactors = FALSE
    )
    colnames(sample_df)[2] <- data_type_label
    filepath <- file.path(outputs_dir, filename)
    save_csv_no_overwrite(sample_df, filepath)
  }
}

message("Exporting per-sample raw counts...")
export_per_sample(raw_counts_df, sample_meta, "raw_counts", outputs_dir)

message("Exporting per-sample normalized counts...")
export_per_sample(norm_counts_df, sample_meta, "DESeq2_normalized", outputs_dir)

message("Exporting per-sample VST values...")
export_per_sample(vst_df, sample_meta, "VST", outputs_dir)

message("Global RNA-seq DESeq2 pipeline complete.")
sessionInfo()
