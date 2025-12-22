# Targeted DESeq2 analysis for an a priori hypothalamus gene set
# Requirements:
# - Working directory contains:
#     quants/                        # per-sample Salmon quant directories
#     tx_id.txt, gene_id.txt         # transcript-to-gene mapping inputs
#     sampleInfo.csv (or sampleInfo_2.csv) with SampleID, condition, sex
#     DeSeq2_hypothalamus_genes.csv # 1-column 'gene' with targeted gene IDs

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("tximport", "DESeq2", "apeglm"), ask = FALSE, update = FALSE)
  install.packages(c("dplyr", "ggplot2", "tidyr", "tibble",
                     "readr", "ggsci", "ggrepel"),
                   quiet = TRUE)
  library(dplyr)
  library(tximport)
  library(DESeq2)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggsci)
  library(ggrepel)
  library(apeglm)
})

# --------------------- Path validation ---------------------
if (!all(file.exists(c("tx_id.txt", "gene_id.txt")))) {
  stop("tx_id.txt or gene_id.txt missing.")
}
if (!dir.exists("quants")) stop("No 'quants/' directory found.")

# --------------------- tx2gene mapping ---------------------
tx_raw   <- read.csv("tx_id.txt", header = FALSE)$V1
gene_raw <- read.csv("gene_id.txt", header = FALSE)$V1
tx2gene  <- data.frame(
  TXNAME = gsub("\\.\\d+$", "", tx_raw),
  GENEID = gene_raw
)
cat("tx2gene ready:", nrow(tx2gene), "transcripts\n")

# --------------------- Import Salmon quantifications ---------------------
quants_dir <- "quants"
quant_dirs <- list.files(quants_dir)
if (length(quant_dirs) == 0) stop("No directories in 'quants/'.")
files <- file.path(quants_dir, quant_dirs, "quant.sf")
if (!all(file.exists(files))) {
  missing <- files[!file.exists(files)]
  stop("Missing quant.sf files: ", paste(head(missing, 3), collapse = ", "))
}
sampleID <- gsub("_quant", "", quant_dirs)
names(files) <- sampleID

cat("Found", length(sampleID), "samples:", paste(head(sampleID, 5), collapse = ", "), "...\n")

# --------------------- tximport ---------------------
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
cat("tximport complete. Counts dimensions:", dim(txi$counts), "\n")

# --------------------- Sample metadata ---------------------
if (file.exists("sampleInfo.csv")) {
  csv_file <- "sampleInfo.csv"
} else if (file.exists("sampleInfo_2.csv")) {
  csv_file <- "sampleInfo_2.csv"
} else {
  stop("No sampleInfo CSV found.")
}
treat <- read.csv(csv_file, header = TRUE)
cat("Using metadata:", csv_file, ". Dimensions:", dim(treat), "\n")

id_col_name <- ifelse("SampleID" %in% names(treat), "SampleID", names(treat)[1])
condition_col_name <- ifelse(
  "condition" %in% names(treat), "condition",
  ifelse("Comments" %in% names(treat), "Comments",
         stop("No 'condition' or 'Comments' column."))
)
sex_col_name <- "sex"

id_col <- treat[[id_col_name]]
cat("Raw CSV IDs (first 5):", paste(head(as.character(id_col), 5), collapse = ", "), "\n")
cat("Quants sampleID (first 5):", paste(head(sampleID, 5), collapse = ", "), "\n")

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

cat("Cleaned CSV IDs (first 5):", paste(head(treat_ids_clean, 5), collapse = ", "), "\n")

condition_col <- treat[[condition_col_name]]
sex_col       <- treat[[sex_col_name]]

samples <- data.frame(
  sampleID  = sampleID,
  condition = condition_col[match(sampleID, treat_ids_clean)],
  sex       = sex_col[match(sampleID, treat_ids_clean)]
)
samples$sex_condition <- factor(paste0(samples$sex, samples$condition))

na_count <- sum(is.na(samples$condition) | is.na(samples$sex))
cat("Successful matches:", length(sampleID) - na_count, "/", length(sampleID), "\n")
print(table(samples$sex_condition, useNA = "ifany"))
if (na_count > 0) {
  mismatches <- sampleID[is.na(samples$condition)]
  stop(
    "Metadata mismatches (", na_count, "): e.g., ",
    paste(head(mismatches, 3), collapse = ", "),
    ". Check ID cleaning logic matches your sampleInfo.csv format."
  )
}
cat("Metadata loaded successfully—proceeding to DESeq2.\n")

# --------------------- Factor conversion ---------------------
samples$condition <- factor(samples$condition, levels = c(0, 1), labels = c("control", "treatment"))
samples$sex <- factor(samples$sex)
cat("Factored metadata:\n")
print(table(samples$sex, samples$condition))

# --------------------- Build DESeqDataSet ---------------------
dds_treat <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ sex + condition)

keep <- rowSums(counts(dds_treat)) >= 10
dds_treat <- dds_treat[keep, ]
cat("Treatment model filtered genes:", dim(dds_treat), "\n")

# --------------------- Load hypothalamus gene list ---------------------
if (!file.exists("DeSeq2_hypothalamus_genes.csv")) {
  stop("DeSeq2_hypothalamus_genes.csv missing.")
}
hypothalamus_genes_csv <- read.csv("DeSeq2_hypothalamus_genes.csv")
keep_hypothalamus_vector <- hypothalamus_genes_csv$gene
keep_hypothalamus_vector <- unique(
  keep_hypothalamus_vector[
    !is.na(keep_hypothalamus_vector) & nzchar(keep_hypothalamus_vector)
  ]
)
cat("Loaded", length(keep_hypothalamus_vector), "hypothalamus genes.\n")

# --------------------- Subset to hypothalamus genes ---------------------
present_hypo_genes <- rownames(dds_treat)[rownames(dds_treat) %in% keep_hypothalamus_vector]
if (length(present_hypo_genes) == 0) {
  cat("Sample dataset genes:", paste(head(rownames(dds_treat), 5), collapse = ", "), "\n")
  cat("Sample hypothalamus genes:", paste(head(keep_hypothalamus_vector, 5), collapse = ", "), "\n")
  stop("No hypothalamus genes found. Check naming (symbols vs IDs).")
}
cat("Found", length(present_hypo_genes), "matching hypothalamus genes in dataset.\n")

dds_hypothalamus <- dds_treat[present_hypo_genes, ]
cat("Hypothalamus subset dimensions:", dim(dds_hypothalamus), "\n")

# --------------------- DESeq on hypothalamus subset ---------------------
set.seed(1)
dds_hypothalamus <- DESeq(dds_hypothalamus, quiet = TRUE)
cat("DESeq on hypothalamus subset complete.\n")

# --------------------- Extract results ---------------------
res_hypothalamus <- results(dds_hypothalamus, alpha = 0.05)
res_hypothalamus_df <- as.data.frame(res_hypothalamus) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)

coefs <- resultsNames(dds_hypothalamus)
coef_name <- if ("condition_treatment_vs_control" %in% coefs) {
  "condition_treatment_vs_control"
} else {
  coefs[grepl("^condition_", coefs)][1]
}
cat("Using coefficient/contrast:", coef_name, "\n")

resLFC_hypothalamus <- lfcShrink(
  dds_hypothalamus,
  res  = res_hypothalamus,
  coef = coef_name,
  type = "apeglm"
)
resLFC_hypothalamus_df <- as.data.frame(resLFC_hypothalamus) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(log2FoldChange)

sig_hypo <- res_hypothalamus_df %>%
  dplyr::filter(!is.na(padj) & padj < 0.05)
cat("Hypothalamus DEGs (padj < 0.05):", nrow(sig_hypo), "\n")
cat("Upregulated:", sum(sig_hypo$log2FoldChange > 0, na.rm = TRUE),
    "; Downregulated:", sum(sig_hypo$log2FoldChange < 0, na.rm = TRUE), "\n")

# --------------------- Save results (timestamped) ---------------------
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

write.csv(
  res_hypothalamus_df,
  paste0("DeSeq2_hypothalamus_results_", ts, ".csv"),
  row.names = FALSE
)
saveRDS(res_hypothalamus, paste0("res_hypothalamus_", ts, ".rds"))

write.csv(
  resLFC_hypothalamus_df,
  paste0("resLFC_hypothalamus_", ts, ".csv"),
  row.names = FALSE
)
saveRDS(resLFC_hypothalamus, paste0("resLFC_hypothalamus_", ts, ".rds"))

saveRDS(dds_hypothalamus, paste0("dds_hypothalamus_", ts, ".rds"))

cat("Results saved with timestamp:", ts, "\n")

# --------------------- Volcano plot ---------------------
volc_hypo <- res_hypothalamus_df %>%
  mutate(
    neg_log_p = -log10(pvalue),
    up_sig    = padj < 0.05 & log2FoldChange > 1,
    down_sig  = padj < 0.05 & log2FoldChange < -1,
    direction = dplyr::case_when(
      up_sig   ~ "Up (sig)",
      down_sig ~ "Down (sig)",
      TRUE     ~ "Not sig"
    ),
    label = dplyr::case_when(
      down_sig & (abs(log2FoldChange) > 1.5 | padj < 0.01) ~ gene,
      TRUE ~ NA_character_
    )
  )

p_volcano <- ggplot(volc_hypo, aes(x = log2FoldChange, y = neg_log_p, color = direction)) +
  geom_point(alpha = ifelse(volc_hypo$direction == "Not sig", 0.4, 0.8), size = 1.8) +
  scale_color_manual(
    values = c("Down (sig)" = "#1f77b4",
               "Up (sig)"   = "#d62728",
               "Not sig"    = "grey60")
  ) +
  geom_text_repel(
    data = subset(volc_hypo, !is.na(label)),
    aes(label = label),
    max.overlaps = Inf,
    size = 3,
    box.padding = 0.4,
    point.padding = 0.4
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.6) +
  labs(
    title = "Hypothalamus Gene Set: Treatment vs Control",
    subtitle = paste0("Coefficient: ", coef_name,
                      " | padj < 0.05 & |log2FC| ≥ 1"),
    x = expression(log[2] ~ "Fold Change"),
    y = expression(-log[10] ~ "(p-value)"),
    color = "Significance"
  ) +
  theme_classic()

ggsave(paste0("volcano_hypothalamus_", ts, ".pdf"),
       p_volcano, width = 7.5, height = 5.5, dpi = 300)
ggsave(paste0("volcano_hypothalamus_", ts, ".png"),
       p_volcano, width = 7.5, height = 5.5, dpi = 600)

# --------------------- Enrichment-ready tables ---------------------
rank_unshrunk <- res_hypothalamus_df %>%
  transmute(
    gene,
    log2FoldChange,
    pvalue,
    padj,
    baseMean
  )
write.csv(
  rank_unshrunk,
  paste0("rank_hypothalamus_unshrunk_for_GSEA_", ts, ".csv"),
  row.names = FALSE
)

rank_shrunk <- resLFC_hypothalamus_df %>%
  transmute(
    gene,
    log2FoldChange = log2FoldChange
  )
write.csv(
  rank_shrunk,
  paste0("rank_hypothalamus_shrunk_for_GSEA_", ts, ".csv"),
  row.names = FALSE
)

cat("Enrichment tables saved. Ready for downstream enrichment.\n")

# --------------------- Supplementary tables (optional) ---------------------
res_unshrunk <- as.data.frame(res_hypothalamus) %>%
  rownames_to_column("gene") %>%
  rename(
    baseMean              = baseMean,
    log2FoldChange_unsh   = log2FoldChange,
    lfcSE_unsh            = lfcSE,
    pvalue                = pvalue,
    padj                  = padj
  )

res_shrunk <- as.data.frame(resLFC_hypothalamus) %>%
  rownames_to_column("gene") %>%
  rename(
    log2FoldChange_shrunk = log2FoldChange,
    lfcSE_shrunk          = lfcSE
  )

supp_table_all <- res_unshrunk %>%
  left_join(res_shrunk %>% select(gene, log2FoldChange_shrunk, lfcSE_shrunk),
            by = "gene") %>%
  transmute(
    Gene                      = gene,
    Base_Mean                 = baseMean,
    Log2FC_Unshrunken         = log2FoldChange_unsh,
    LFC_SE_Unshrunken         = lfcSE_unsh,
    Unadjusted_p              = pvalue,
    Adjusted_p_FDR            = padj,
    Log2FC_Shrunken_apeglm    = log2FoldChange_shrunk,
    LFC_SE_Shrunken_apeglm    = lfcSE_shrunk
  ) %>%
  arrange(Adjusted_p_FDR, Unadjusted_p)

ts_full <- format(Sys.time(), "%Y%m%d")
write.csv(
  supp_table_all,
  paste0("Supplement_Table_Hypothalamus_All_Genes_", ts_full, ".csv"),
  row.names = FALSE,
  na = ""
)
readr::write_tsv(
  supp_table_all,
  paste0("Supplement_Table_Hypothalamus_All_Genes_", ts_full, ".tsv")
)

cat("Hypothalamus analysis complete.\n")
sessionInfo()
