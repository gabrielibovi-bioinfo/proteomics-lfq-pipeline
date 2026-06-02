################################################################################
# Script: 02_dep_analysis.R
# Description: Differential expression analysis visualization for LFQ
#              proteomics data. Generates volcano plot, heatmap of DEPs,
#              and PCA from limma results.
#
# Steps performed:
#   1. Load limma results (MO_vs_SC contrast)
#   2. Suggest logFC thresholds based on percentiles
#   3. logFC density plot
#   4. Categorize proteins as up, down, or not significant
#   5. Volcano plot with top 5 labeled proteins
#   6. Heatmap of differentially expressed proteins (DEPs)
#   7. PCA of DEPs
#
# Input files:
#   - limma_rlr_BH/per_contrast_results/MO_vs_SC.csv
#
# Output files:
#   - analise_DEP/logFC_percentil.xlsx
#   - analise_DEP/genes_lFC.xlsx
#   - analise_DEP/contagem_degs_sig.xlsx
#   - analise_DEP/01_volcano_plot.tiff
#   - analise_DEP/02_heatmap_DEP.tiff
#   - analise_DEP/03_pca_pvalue.tiff
################################################################################

library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(readr)
library(tidyr)
library(tibble)
library(readxl)
library(ggrepel)
library(gprofiler2)
library(clusterProfiler)
library(enrichplot)
library(writexl)

options(scipen = 999)
set.seed(123)

# ==============================================================================
# SECTION 1: Load limma results
# ==============================================================================

dados <- read_csv("limma_rlr_BH/per_contrast_results/MO_vs_SC.csv", show_col_types = FALSE)

summary(dados)
summary(dados$logFC)

# ==============================================================================
# SECTION 2: Suggest logFC thresholds based on percentiles
# ==============================================================================

suggest_logfc_thresholds <- function(logfc_vector,
                                     percentiles = c(0.2, 0.4, 0.6, 0.8, 1.0)) {
  logfc_abs  <- abs(logfc_vector[!is.na(logfc_vector)])
  thresholds <- quantile(logfc_abs, probs = percentiles)
  counts     <- sapply(thresholds, function(th) sum(logfc_abs > th))

  data.frame(
    Percentile       = percentiles * 100,
    logFC_Threshold  = as.numeric(thresholds),
    Num_Proteins     = counts
  )
}

logfc_thresholds <- suggest_logfc_thresholds(dados$logFC)
print(logfc_thresholds)
# write_xlsx(logfc_thresholds, "analise_DEP/logFC_percentil.xlsx")

# NOTE: A logFC percentile of 90 means 90% of proteins have |logFC| below that
# threshold and 10% above. Use this to calibrate your significance cutoff.

# ==============================================================================
# SECTION 3: logFC density plot
# ==============================================================================

lfc_density <- ggplot(dados, aes(x = logFC)) +
  geom_density(fill = "#A1C9F4", alpha = 0.6) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  geom_vline(xintercept = c(-0.25, 0.25), color = "red", linetype = "dotted") +
  labs(
    title = "logFC Density",
    x     = "log Fold Change",
    y     = "Density"
  ) +
  theme_minimal(base_size = 13)

print(lfc_density)
# ggsave(plot = lfc_density, "analise_DEP/LFC_density.tiff", width = 5, height = 4, dpi = 300, bg = "white")

# ==============================================================================
# SECTION 4: Categorize proteins by expression direction
# ==============================================================================

# logFC threshold set to 0 due to small variance between treatment groups
logfc_threshold <- 0
pval_threshold  <- 0.05
                       
dados <- dados %>%
  mutate(
    expression_category = case_when(
      P.Value < pval_threshold & logFC < 0 ~ "down",
      P.Value < pval_threshold & logFC > 0 ~ "up",
      TRUE ~ "NS"
    )
  )

# write_xlsx(dados, "analise_DEP/genes_lFC.xlsx")

dep_counts <- dados %>%
  count(expression_category) %>%
  arrange(desc(n))

print(dep_counts)
# write_xlsx(dep_counts, "analise_DEP/counts_degs_sig.xlsx")

# ==============================================================================
# SECTION 5: Volcano plot
# ==============================================================================

top5_proteins <- dados %>%
  filter(P.Value < pval_threshold, abs(logFC) > logfc_threshold) %>%
  arrange(P.Value, desc(abs(logFC))) %>%
  head(5)

volcano_plot <- ggplot(dados,
                       aes(x = logFC, y = -log10(P.Value),
                           color = expression_category)) +
  geom_point(alpha = 0.7, shape = 19, size = 2.4) +
  scale_color_manual(values = c("up" = "olivedrab", "down" = "darkred", "NS" = "black")) +
  geom_text_repel(
    data        = top5_proteins,
    aes(label   = uniprot_id),
    size        = 4,
    box.padding = 0.5,
    max.overlaps = 10
  ) +
  geom_hline(yintercept = -log10(pval_threshold),
             linetype   = "dashed", color = "gray40") +
  theme_minimal() +
  labs(
    title = "Volcano Plot — MO vs SC",
    x     = "Log Fold Change",
    y     = "-Log10 P-Value",
    color = "Expression"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

print(volcano_plot)
# ggsave(plot = volcano_plot, "analise_DEP/01_volcano_plot.tiff", width = 6, height = 5, dpi = 300, bg = "white")

# ==============================================================================
# SECTION 6: Heatmap of differentially expressed proteins (DEPs)
# ==============================================================================

# Filter to DEPs only; select LFQ intensity columns
deps <- dados %>%
  filter(P.Value < pval_threshold & abs(logFC) > logfc_threshold) %>%
  dplyr::select(uniprot_id, starts_with("LFQ intensity")) %>%
  column_to_rownames(var = "uniprot_id")

# NOTE: Change P.Value to adj.P.Val above for a stricter threshold

lfq_cols_dep     <- grep("^LFQ intensity", colnames(deps))
dados_lfq_dep    <- deps[, lfq_cols_dep]

sample_groups    <- ifelse(grepl("MO", colnames(dados_lfq_dep)), "MO", "SC")
annotation_col   <- data.frame(Group = sample_groups)
rownames(annotation_col) <- colnames(dados_lfq_dep)

group_colors_dep <- list(Group = c("MO" = "seagreen", "SC" = "chocolate"))
heatmap_colors   <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(100)

dep_matrix <- as.matrix(dados_lfq_dep)

heatmap_dep <- pheatmap(
  mat                      = dep_matrix,
  scale                    = "row",      # Options: "row", "column", "none"
  border_color             = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D",
  annotation_col           = annotation_col,
  annotation_colors        = group_colors_dep,
  color                    = heatmap_colors,
  main                     = "LFQ Intensity — DEPs (row-scaled)",
  show_rownames            = TRUE,
  cluster_cols             = TRUE,
  fontsize_row             = 8,
  fontsize_col             = 8,
  na_col                   = "#053061"
)

# tiff("analise_DEP/02_heatmap_DEP.tiff", width = 7, height = 7, units = "in", res = 300, bg = "white")
# print(heatmap_dep)
# dev.off()

# ==============================================================================
# SECTION 7: PCA of differentially expressed proteins
# ==============================================================================

dados_log <- log2(dados_lfq_dep + 1)
dados_log <- dados_log[rowSums(is.na(dados_log)) == 0, ]
dados_log <- dados_log[rowMeans(is.na(dados_log)) < 0.5, ]

pca_res <- prcomp(t(dados_log), scale. = TRUE)

pca_df        <- as.data.frame(pca_res$x)
pca_df$Group  <- ifelse(grepl("MO", rownames(pca_df)), "MO", "SC")

pca_dep <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 5) +
  theme_bw() +
  labs(
    title = "PCA — LFQ Intensity (DEPs, P-value)",
    x     = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
    y     = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)")
  ) +
  scale_color_manual(values = c("MO" = "seagreen", "SC" = "chocolate"))

print(pca_dep)
# ggsave(plot = pca_dep, "analise_DEP/03_pca_pvalue.tiff", width = 5.5, height = 4.4, dpi = 300, bg = "white")
