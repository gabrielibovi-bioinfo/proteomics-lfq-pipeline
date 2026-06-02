################################################################################
# Script: 01_qc_normalization.R
# Description: Quality control, preprocessing, normalization, and differential
#              expression analysis of LFQ proteomics data using proteoDA.
#
# Steps performed:
#   1. Load raw LFQ intensity data from MaxQuant output (Excel)
#   2. Preprocess: log2 transform, remove zero/constant proteins
#   3. Quality control: heatmaps, Pearson correlation matrix, PCA (all samples)
#   4. Remove outlier sample (SC-1) and repeat QC
#   5. Normalize with proteoDA (RLR method) and generate normalization report
#   6. Fit limma model and extract differentially expressed proteins (DEPs)
#   7. Export limma results and diagnostic plots
#
# Input files:
#   - data/proteomics_data.xlsx    (MaxQuant output with LFQ intensity columns)
#
# Output files:
#   - Normalization_report/        (normalization and QC PDF reports)
#   - limma_rlr_BH/                (limma results per contrast)
#   - figures/pre-processing/      (heatmaps, PCA plots)
#
# NOTE: Update FILE_PATH at the top of this script before running.
################################################################################

library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(pheatmap)
library(factoextra)
library(plotly)
library(ggfortify)
library(proteoDA)

options(scipen = 999)
set.seed(123)

# NOTE: Install dependencies (run once):
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("preprocessCore", "affy", "vsn"))
# devtools::install_github("ByrumLab/proteoDA", dependencies = TRUE,
#                          build_vignettes = TRUE, force = TRUE)

# ==============================================================================
# SECTION 0: File path — update before running
# ==============================================================================

FILE_PATH <- "data/proteomics_data.xlsx"   # Update to your file path

# ==============================================================================
# SECTION 1: Load data and extract LFQ intensity columns
# ==============================================================================

dados    <- read_excel(FILE_PATH)
lfq_cols <- grep("^LFQ intensity", colnames(dados))
dados_lfq <- dados[, lfq_cols]

# ==============================================================================
# SECTION 2: Preprocessing — log2 transform, filter, remove constant proteins
# ==============================================================================

# Log2 transform (add 1 to avoid log(0))
dados_lfq <- log2(dados_lfq + 1)

# Remove proteins with all-zero values
dados_lfq <- dados_lfq[rowSums(dados_lfq > 0) > 0, ]

# Remove proteins with all missing values
dados_lfq <- dados_lfq[rowSums(!is.na(dados_lfq)) >= 1, ]

# Remove constant proteins (zero variance)
vars <- apply(dados_lfq, 1, var, na.rm = TRUE)
dados_lfq <- dados_lfq[vars > 0, ]

cat("Missing values after preprocessing:", sum(is.na(dados_lfq)), "\n")
cat("Proteins remaining:", nrow(dados_lfq), "\n")

# Group annotation for heatmap/PCA
sample_groups  <- ifelse(grepl("MO", colnames(dados_lfq)), "MO", "SC")
annotation_col <- data.frame(Group = sample_groups)
rownames(annotation_col) <- colnames(dados_lfq)
group_colors   <- list(Group = c("MO" = "chocolate", "SC" = "seagreen"))

# ==============================================================================
# SECTION 3: Quality control — all samples (before outlier removal)
# ==============================================================================

# --- 3.1: Heatmap — all proteins ---
heatmap_all <- pheatmap(
  mat                      = dados_lfq,
  scale                    = "row",
  border_color             = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D",
  annotation_col           = annotation_col,
  annotation_colors        = group_colors,
  color                    = colorRampPalette(c("darkblue", "white", "darkred"))(50),
  main                     = "LFQ Intensity — All samples",
  show_rownames            = FALSE,
  cluster_cols             = FALSE
)
# ggsave(plot = heatmap_all, "figures/pre-processing/01_heatmap_all.tiff", width = 5.1, height = 5.3, dpi = 300, bg = "white")

# --- 3.2: Pearson correlation heatmap ---
cor_matrix    <- cor(dados_lfq, use = "pairwise.complete.obs", method = "pearson")

heatmap_corr <- pheatmap(
  mat                      = cor_matrix,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D",
  annotation_col           = annotation_col,
  annotation_colors        = group_colors,
  color                    = colorRampPalette(c("darkblue", "white", "darkred"))(25),
  main                     = "Pearson Correlation — LFQ Intensity",
  show_rownames            = TRUE,
  show_colnames            = FALSE,
  border_color             = FALSE
)
# ggsave(plot = heatmap_corr, "figures/pre-processing/02_pearson_correlation.tiff", width = 5, height = 5, dpi = 300, bg = "white")

# --- 3.3: PCA — all samples ---
pca_res     <- prcomp(t(dados_lfq), center = TRUE, scale. = TRUE)
var_exp     <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
pc1_pct     <- round(var_exp[1], 2)
pc2_pct     <- round(var_exp[2], 2)

pca_df      <- as.data.frame(pca_res$x)
pca_df$Sample <- rownames(pca_res$x)
pca_df$Group  <- ifelse(grepl("MO", rownames(pca_res$x)), "MO", "SC")

pca_all <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(aes(fill = Group), color = "white", shape = 21, size = 8, stroke = 1) +
  theme_minimal() +
  labs(
    title = "PCA — LFQ Intensity (all samples)",
    x     = paste0("PC1 (", pc1_pct, "%)"),
    y     = paste0("PC2 (", pc2_pct, "%)")
  ) +
  scale_color_manual(values = c("MO" = "chocolate", "SC" = "seagreen")) +
  theme(
    axis.title  = element_text(size = 14),
    axis.text   = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12)
  )

print(pca_all)
# ggsave(plot = pca_all, "figures/pre-processing/03_pca_all.tiff", width = 5.5, height = 4.4, dpi = 300, bg = "white")

# ==============================================================================
# SECTION 4: Remove outlier sample (SC-1) and repeat QC
# ==============================================================================

dados_lfq_filt  <- dados_lfq %>% select(-`LFQ intensity SC-1`)
dados_lfq_filt  <- dados_lfq_filt[rowSums(dados_lfq_filt > 0) > 0, ]

sample_groups_f <- ifelse(grepl("MO", colnames(dados_lfq_filt)), "MO", "SC")
annot_filt      <- data.frame(Group = sample_groups_f)
rownames(annot_filt) <- colnames(dados_lfq_filt)

# --- 4.1: Heatmap — after outlier removal ---
heatmap_filt <- pheatmap(
  mat                      = dados_lfq_filt,
  scale                    = "row",
  border_color             = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D",
  annotation_col           = annot_filt,
  annotation_colors        = group_colors,
  color                    = colorRampPalette(c("darkblue", "white", "darkred"))(25),
  main                     = "LFQ Intensity — After outlier removal",
  show_rownames            = FALSE,
  cluster_cols             = FALSE
)
# ggsave(plot = heatmap_filt, "figures/pre-processing/04_heatmap_no_SC1.tiff", width = 5.1, height = 5.3, dpi = 300, bg = "white")

# --- 4.2: Pearson correlation — after outlier removal ---
cor_filt      <- cor(dados_lfq_filt, use = "pairwise.complete.obs", method = "pearson")

heatmap_filt_corr <- pheatmap(
  mat                      = cor_filt,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D",
  annotation_col           = annot_filt,
  annotation_colors        = group_colors,
  color                    = colorRampPalette(c("darkblue", "white", "darkred"))(25),
  main                     = "Pearson Correlation — After outlier removal",
  show_rownames            = TRUE,
  show_colnames            = FALSE,
  border_color             = FALSE
)
# ggsave(plot = heatmap_filt_corr, "figures/pre-processing/05_pearson_no_SC1.tiff", width = 5, height = 5, dpi = 300, bg = "white")

# --- 4.3: PCA — after outlier removal ---
pca_filt_res  <- prcomp(t(dados_lfq_filt), center = TRUE, scale. = TRUE)
var_exp_f     <- (pca_filt_res$sdev^2) / sum(pca_filt_res$sdev^2) * 100

pca_filt_df   <- as.data.frame(pca_filt_res$x)
pca_filt_df$Sample <- rownames(pca_filt_res$x)
pca_filt_df$Group  <- ifelse(grepl("MO", rownames(pca_filt_res$x)), "MO", "SC")

pca_filt <- ggplot(pca_filt_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(aes(fill = Group), color = "white", shape = 21, size = 8, stroke = 1) +
  theme_minimal() +
  labs(
    title = "PCA — LFQ Intensity (after outlier removal)",
    x     = paste0("PC1 (", round(var_exp_f[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_exp_f[2], 2), "%)")
  ) +
  scale_color_manual(values = c("MO" = "chocolate", "SC" = "seagreen")) +
  theme(
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12)
  )

print(pca_filt)
# ggsave(plot = pca_filt, "figures/pre-processing/06_pca_no_SC1.tiff", width = 5.5, height = 4.4, dpi = 300, bg = "white")

# ==============================================================================
# SECTION 5: Normalization with proteoDA (RLR)
# ==============================================================================

# Reload full dataset for proteoDA workflow
dados_raw  <- read_excel(FILE_PATH)
lfq_cols   <- grep("^LFQ intensity", colnames(dados_raw))
colnames(dados_raw)[colnames(dados_raw) == "T: Gene names"] <- "uniprot_id"
dados_raw  <- dados_raw[!is.na(dados_raw$uniprot_id), ]

dados_lfq2 <- as.data.frame(dados_raw[, lfq_cols])
rownames(dados_lfq2) <- dados_raw$uniprot_id

# Remove outlier sample
dados_lfq2 <- dados_lfq2[, colnames(dados_lfq2) != "LFQ intensity SC-1"]

# Filter proteins: remove all-zero rows and constant rows
dados_lfq2 <- dados_lfq2[rowSums(dados_lfq2 > 0) > 0, ]
dados_lfq2 <- dados_lfq2[rowSums(!is.na(dados_lfq2)) >= 1, ]
vars2      <- apply(dados_lfq2, 1, var, na.rm = TRUE)
dados_lfq2 <- dados_lfq2[vars2 > 0, ]

# Build sample metadata
sample_groups2   <- ifelse(grepl("MO", colnames(dados_lfq2)), "MO", "SC")
sample_metadata  <- data.frame(
  data_column_name = colnames(dados_lfq2),
  group            = sample_groups2
)
rownames(sample_metadata) <- colnames(dados_lfq2)
sample_metadata <- sample_metadata[sample_metadata$data_column_name != "LFQ intensity SC-1", ]

# Build annotation table
annotation <- data.frame(uniprot_id = rownames(dados_lfq2))

# Build DAList object
raw <- DAList(
  data       = as.matrix(dados_lfq2),
  annotation = annotation,
  metadata   = sample_metadata
)

# Filter proteins present in >= 50% of samples in at least one group
filtered <- raw |>
  zero_to_missing() |>
  filter_proteins_by_proportion(min_prop = 0.5, grouping_column = "group")

# Generate normalization report
write_norm_report(
  filtered,
  grouping_column      = "group",
  output_dir           = "Normalization_report",
  filename             = NULL,
  overwrite            = TRUE,
  suppress_zoom_legend = FALSE,
  use_ggrastr          = FALSE
)

# Normalize using RLR (robust linear regression)
normalized <- normalize_data(filtered, norm_method = "rlr") # Other options: "cycloess", "vsn", "quantile"

# Generate QC report
write_qc_report(
  normalized,
  color_column      = "group",
  label_column      = NULL,
  output_dir        = "Normalization_report",
  filename          = "QC_Report_rlr.pdf",
  overwrite         = TRUE,
  top_proteins      = 500,
  standardize       = TRUE,
  pca_axes          = c(1, 2),
  dist_metric       = "euclidean",
  clust_method      = "complete",
  show_all_proteins = FALSE
)

# Normalization method comparison (number of DEPs found):
#   cycloess | BH: 11 up - 13 down = 24
#   vsn      | BH: 11 up - 12 down = 23
#   rlr      | BH: 11 up - 15 down = 26  *** Best overall performance
#   quantile | BH:  9 up - 14 down = 23

# ==============================================================================
# SECTION 6: Differential expression analysis with limma
# ==============================================================================

normalized$metadata$group <- factor(normalized$metadata$group, levels = c("MO", "SC"))

# Build model without intercept
no_intercept <- add_design(normalized, design_formula = ~0 + group)

# Define contrast: MO vs SC
no_intercept <- add_contrasts(no_intercept, contrasts_vector = c("MO_vs_SC = MO - SC"))

# Fit limma model
fit <- fit_limma_model(no_intercept)

# Extract DEP results
results <- extract_DA_results(
  fit,
  pval_thresh       = 0.05,
  lfc_thresh        = 0.25,
  adj_method        = "BH",       # Options: "BH", "BY", "holm", "none"
  extract_intercept = FALSE
)

# ==============================================================================
# SECTION 7: Export results
# ==============================================================================

write_limma_tables(
  results,
  output_dir        = "limma_rlr_BH",
  overwrite         = TRUE,
  contrasts_subdir  = "per_contrast_results",
  summary_csv       = NULL,
  combined_file_csv = NULL,
  spreadsheet_xlsx  = NULL,
  add_filter        = TRUE
)

write_limma_plots(
  results,
  output_dir       = "limma_rlr_BH",
  grouping_column  = "group",
  overwrite        = TRUE
)
