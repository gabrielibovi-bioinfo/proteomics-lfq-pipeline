################################################################################
# Script: 03_gsea.R
# Description: Gene Set Enrichment Analysis (GSEA) for LFQ proteomics data
#              using GO, KEGG, and Reactome databases.
#              Organism: Rattus norvegicus (rat)
#              Contrast: MO vs SC
#
# Steps performed:
#   1. Load limma results and map gene symbols to Entrez IDs
#   2. Build ranked gene list using limma t-statistic
#   3. GSEA — Gene Ontology (GO): BP, MF, CC
#   4. GSEA — KEGG pathways
#   5. GSEA — Reactome pathways
#   Each section generates: dotplot, heatplot, ridgeplot, GSEA plots,
#   enrichment map, and network plot (cnetplot)
#
# Input files:
#   - limma_rlr_BH/per_contrast_results/MO_vs_SC.csv
#
# Output files:
#   - gsea_GO/          (GO results and figures)
#   - gsea_kegg/        (KEGG results and figures)
#   - gsea_reactome/    (Reactome results and figures)
################################################################################

library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(ReactomePA)
library(dplyr)
library(ggplot2)
library(writexl)

# NOTE: Install dependencies (run once):
# BiocManager::install(c("clusterProfiler", "org.Rn.eg.db", "enrichplot",
#                         "AnnotationDbi", "ReactomePA"))

options(scipen = 999)
set.seed(123)

# ==============================================================================
# SECTION 1: Load limma results and build ranked gene list
# ==============================================================================

limma_res <- read.csv("limma_rlr_BH/per_contrast_results/MO_vs_SC.csv",
                      stringsAsFactors = FALSE)

# Filter missing values and remove duplicates
limma_res <- limma_res %>%
  filter(!is.na(uniprot_id), !is.na(logFC)) %>%
  distinct(uniprot_id, .keep_all = TRUE)

# Map gene symbols to Entrez IDs (rat genome)
gene_map  <- bitr(limma_res$uniprot_id,
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Rn.eg.db)

limma_res <- merge(limma_res, gene_map,
                   by.x = "uniprot_id", by.y = "SYMBOL")

# Build ranked list using t-statistic (preferred over logFC for GSEA)
gene_list         <- limma_res$t
names(gene_list)  <- limma_res$ENTREZID
gene_list         <- na.omit(gene_list)
gene_list         <- sort(gene_list, decreasing = TRUE)

cat("Genes in ranked list:", length(gene_list), "\n")

# ==============================================================================
# SECTION 2: GSEA — Gene Ontology (GO)
# ==============================================================================

gsea_go <- gseGO(
  geneList       = gene_list,
  OrgDb          = org.Rn.eg.db,
  ont            = "ALL",       # BP (Biological Process), MF, CC, or ALL
  keyType        = "ENTREZID",
  minGSSize      = 5,
  maxGSSize      = 800,
  pvalueCutoff   = 0.05,
  pAdjustMethod  = "none",      # Options: "none", "BH", "bonferroni"
  verbose        = TRUE
)

gsea_go <- setReadable(gsea_go, OrgDb = org.Rn.eg.db)

# Export GO results
gsea_go_df <- gsea_go@result
# write.csv(gsea_go_df, "gsea_GO/GSEA_GO_MO_vs_SC.csv", row.names = FALSE)

# --- 2.1: Dotplot — faceted by ontology ---
go_dot1 <- dotplot(
  gsea_go, showCategory = 10, split = "ONTOLOGY",
  title = "Enriched Pathways — GO", color = "pvalue",
  font.size = 12, label_format = 50
) +
  facet_grid(. ~ ONTOLOGY, drop = TRUE, scales = "free_x") +
  theme(axis.text.y = element_text(size = 10))

print(go_dot1)
# ggsave("gsea_GO/01_GSEA_GO_ontology.tiff",
#        plot = go_dot1, width = 9, height = 7, dpi = 400, compression = "lzw")

# --- 2.1b: Dotplot — faceted by direction and ontology ---
go_dot2 <- dotplot(
  gsea_go, showCategory = 10, split = ".sign",
  title = "Enriched Pathways — GO (by direction)", color = "pvalue",
  font.size = 12, label_format = 50
) +
  facet_grid(rows = vars(.sign), cols = vars(ONTOLOGY),
             scales = "free_y", drop = TRUE) +
  theme(axis.text.y = element_text(size = 10))

print(go_dot2)
# ggsave("gsea_GO/02_GSEA_GO_direction.tiff",
#        plot = go_dot2, width = 9, height = 6, dpi = 400, compression = "lzw")

# --- 2.2: Heatplot — pathways vs proteins ---
go_heat <- heatplot(gsea_go, foldChange = gene_list, showCategory = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(go_heat)
# ggsave("gsea_GO/03_GSEA_GO_heatplot.tiff",
#        plot = go_heat, width = 10, height = 5, dpi = 400, compression = "lzw")

# --- 2.3: Ridgeplot — enrichment distribution ---
go_ridge <- ridgeplot(gsea_go, showCategory = 15) +
  labs(x = "Enrichment Distribution") +
  theme(
    axis.text    = element_text(size = 12),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

print(go_ridge)
# ggsave("gsea_GO/04_GSEA_GO_ridgeplot.tiff",
#        plot = go_ridge, width = 8, height = 6, dpi = 400, compression = "lzw")

# --- 2.4: GSEA running score plots — top 3 pathways ---
go_gs1 <- gseaplot(gsea_go, by = "all", title = gsea_go$Description[1], geneSetID = 1)
go_gs2 <- gseaplot(gsea_go, by = "all", title = gsea_go$Description[2], geneSetID = 2)
go_gs3 <- gseaplot(gsea_go, by = "all", title = gsea_go$Description[3], geneSetID = 3)
# ggsave("gsea_GO/05_GSEA_GO_top1.tiff", plot = go_gs1, width = 6, height = 5, dpi = 400, compression = "lzw")
# ggsave("gsea_GO/06_GSEA_GO_top2.tiff", plot = go_gs2, width = 6, height = 5, dpi = 400, compression = "lzw")
# ggsave("gsea_GO/07_GSEA_GO_top3.tiff", plot = go_gs3, width = 6, height = 5, dpi = 400, compression = "lzw")

go_gsea1 <- gseaplot2(gsea_go, title = gsea_go$Description[1], geneSetID = 1)
go_gsea2 <- gseaplot2(gsea_go, title = gsea_go$Description[2], geneSetID = 2)
go_gsea3 <- gseaplot2(gsea_go, title = gsea_go$Description[3], geneSetID = 3)
# ggsave("gsea_GO/08_GSEA_GO_top1_v2.tiff", plot = go_gsea1, width = 6, height = 5, dpi = 400, compression = "lzw")
# ggsave("gsea_GO/09_GSEA_GO_top2_v2.tiff", plot = go_gsea2, width = 6, height = 5, dpi = 400, compression = "lzw")
# ggsave("gsea_GO/10_GSEA_GO_top3_v2.tiff", plot = go_gsea3, width = 6, height = 5, dpi = 400, compression = "lzw")

go_gsea5 <- gseaplot2(gsea_go, geneSetID = 1:5,
                      color = c("yellow3", "darkblue", "darkolivegreen", "chocolate", "seagreen"))
# ggsave("gsea_GO/11_GSEA_GO_top5_rank.tiff",
#        plot = go_gsea5, width = 7, height = 5, dpi = 400, compression = "lzw")

# --- 2.5: Network plot (cnetplot) ---
# Shows linkages between genes and enriched GO terms
go_cnet <- cnetplot(gsea_go, showCategory = 10, categorySize = "pvalue",
                    color.params = list(foldChange = gene_list))
# ggsave("gsea_GO/12_GSEA_GO_cnetplot.tiff",
#        plot = go_cnet, width = 12, height = 10, dpi = 400, compression = "lzw")

# --- 2.6: Enrichment map ---
# Overlapping gene sets cluster together, revealing functional modules
gsea_go_sim <- pairwise_termsim(gsea_go)
go_emap     <- emapplot(gsea_go_sim, showCategory = 15)
# ggsave("gsea_GO/13_GSEA_GO_emapplot.tiff",
#        plot = go_emap, width = 10, height = 7, dpi = 400, compression = "lzw")

# ==============================================================================
# SECTION 3: GSEA — KEGG pathways
# ==============================================================================

gsea_kegg <- gseKEGG(
  geneList      = gene_list,
  organism      = "rno",       # rno = Rattus norvegicus
  minGSSize     = 5,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "none",
  verbose       = TRUE
)

gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")

# Export KEGG results
# write.csv(gsea_kegg@result, "gsea_kegg/GSEA_KEGG_MO_vs_SC.csv", row.names = FALSE)

# --- 3.1: Dotplot ---
kegg_dot1 <- dotplot(
  gsea_kegg, showCategory = 10,
  title = "Enriched Pathways — KEGG", color = "pvalue",
  font.size = 12, label_format = 50
) +
  facet_grid(drop = TRUE, scales = "free_x") +
  theme(axis.text.y = element_text(size = 10))
# ggsave("gsea_kegg/01_GSEA_KEGG.tiff", plot = kegg_dot1, width = 7, height = 7, dpi = 400, compression = "lzw")

kegg_dot2 <- dotplot(
  gsea_kegg, showCategory = 10, split = ".sign",
  title = "Enriched Pathways — KEGG (by direction)", color = "pvalue",
  font.size = 12, label_format = 50
) +
  facet_grid(rows = vars(.sign), scales = "free_y", drop = TRUE) +
  theme(axis.text.y = element_text(size = 10))
# ggsave("gsea_kegg/02_GSEA_KEGG_direction.tiff", plot = kegg_dot2, width = 7, height = 6, dpi = 400, compression = "lzw")

# --- 3.2: Heatplot ---
kegg_heat <- heatplot(gsea_kegg, foldChange = gene_list, showCategory = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
# ggsave("gsea_kegg/03_GSEA_KEGG_heatplot.tiff", plot = kegg_heat, width = 10, height = 5, dpi = 400, compression = "lzw")

# --- 3.3: Ridgeplot ---
kegg_ridge <- ridgeplot(gsea_kegg, showCategory = 10) +
  labs(x = "Enrichment Distribution KEGG") +
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
# ggsave("gsea_kegg/04_GSEA_KEGG_ridgeplot.tiff", plot = kegg_ridge, width = 8, height = 6, dpi = 400, compression = "lzw")

# --- 3.4: GSEA running score plots ---
kegg_gs1 <- gseaplot(gsea_kegg, by = "all", title = gsea_kegg$Description[1], geneSetID = 1)
kegg_gs2 <- gseaplot(gsea_kegg, by = "all", title = gsea_kegg$Description[2], geneSetID = 2)
kegg_gs3 <- gseaplot(gsea_kegg, by = "all", title = gsea_kegg$Description[3], geneSetID = 3)
# ggsave("gsea_kegg/05_GSEA_KEGG_top1.tiff", plot = kegg_gs1, width = 6, height = 5, dpi = 400, compression = "lzw")
# ggsave("gsea_kegg/06_GSEA_KEGG_top2.tiff", plot = kegg_gs2, width = 6, height = 5, dpi = 400, compression = "lzw")
# ggsave("gsea_kegg/07_GSEA_KEGG_top3.tiff", plot = kegg_gs3, width = 6, height = 5, dpi = 400, compression = "lzw")

kegg_gsea1 <- gseaplot2(gsea_kegg, title = gsea_kegg$Description[1], geneSetID = 1)
kegg_gsea2 <- gseaplot2(gsea_kegg, title = gsea_kegg$Description[2], geneSetID = 2)
kegg_gsea3 <- gseaplot2(gsea_kegg, title = gsea_kegg$Description[3], geneSetID = 3)
# ggsave("gsea_kegg/08_GSEA_KEGG_top1_v2.tiff", plot = kegg_gsea1, width = 6, height = 5, dpi = 400, compression = "lzw")
# ggsave("gsea_kegg/09_GSEA_KEGG_top2_v2.tiff", plot = kegg_gsea2, width = 6, height = 5, dpi = 400, compression = "lzw")
# ggsave("gsea_kegg/10_GSEA_KEGG_top3_v2.tiff", plot = kegg_gsea3, width = 6, height = 5, dpi = 400, compression = "lzw")

kegg_gsea5 <- gseaplot2(gsea_kegg, geneSetID = 1:5,
                        color = c("yellow3", "darkblue", "darkolivegreen", "chocolate", "seagreen"))
# ggsave("gsea_kegg/11_GSEA_KEGG_top5_rank.tiff", plot = kegg_gsea5, width = 7, height = 5, dpi = 400, compression = "lzw")

# --- 3.5: Network plot ---
kegg_cnet <- cnetplot(gsea_kegg, showCategory = 5, categorySize = "pvalue",
                      color.params = list(foldChange = gene_list))
# ggsave("gsea_kegg/12_GSEA_KEGG_cnetplot.tiff", plot = kegg_cnet, width = 12, height = 10, dpi = 400, compression = "lzw")

# --- 3.6: Enrichment map ---
gsea_kegg_sim <- pairwise_termsim(gsea_kegg)
kegg_emap     <- emapplot(gsea_kegg_sim, showCategory = 10)
# ggsave("gsea_kegg/13_GSEA_KEGG_emapplot.tiff", plot = kegg_emap, width = 10, height = 7, dpi = 400, compression = "lzw")

# ==============================================================================
# SECTION 4: GSEA — Reactome pathways
# ==============================================================================

gsea_reactome <- gsePathway(
  geneList      = gene_list,
  organism      = "rat",
  minGSSize     = 5,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "none",
  verbose       = TRUE
)

gsea_reactome <- setReadable(gsea_reactome, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")

# Export Reactome results
# write.csv(gsea_reactome@result, "gsea_reactome/GSEA_Reactome_MO_vs_SC.csv", row.names = FALSE)

# --- 4.1: Dotplot ---
reac_dot1 <- dotplot(
  gsea_reactome, showCategory = 10,
  title = "Enriched Pathways — Reactome", color = "pvalue",
  font.size = 12, label_format = 50
) +
  facet_grid(drop = TRUE, scales = "free_x") +
  theme(axis.text.y = element_text(size = 10))
# ggsave("gsea_reactome/01_GSEA_Reactome.tiff", plot = reac_dot1, width = 7, height = 7, dpi = 400, compression = "lzw")

reac_dot2 <- dotplot(
  gsea_reactome, showCategory = 10, split = ".sign",
  title = "Enriched Pathways — Reactome (by direction)", color = "pvalue",
  font.size = 12, label_format = 50
) +
  facet_grid(rows = vars(.sign), scales = "free_y", drop = TRUE) +
  theme(axis.text.y = element_text(size = 10))
# ggsave("gsea_reactome/02_GSEA_Reactome_direction.tiff", plot = reac_dot2, width = 7, height = 6, dpi = 400, compression = "lzw")

# --- 4.2: Heatplot ---
reac_heat <- heatplot(gsea_reactome, foldChange = gene_list, showCategory = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
# ggsave("gsea_reactome/03_GSEA_Reactome_heatplot.tiff", plot = reac_heat, width = 10, height = 5, dpi = 400, compression = "lzw")

# --- 4.3: Ridgeplot ---
reac_ridge <- ridgeplot(gsea_reactome, showCategory = 10) +
  labs(x = "Enrichment Distribution Reactome") +
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
# ggsave("gsea_reactome/04_GSEA_Reactome_ridgeplot.tiff", plot = reac_ridge, width = 8, height = 6, dpi = 400, compression = "lzw")

# --- 4.4: GSEA running score plots ---
reac_gs1 <- gseaplot(gsea_reactome, by = "all", title = gsea_reactome$Description[1], geneSetID = 1)
reac_gs2 <- gseaplot(gsea_reactome, by = "all", title = gsea_reactome$Description[2], geneSetID = 2)
reac_gs3 <- gseaplot(gsea_reactome, by = "all", title = gsea_reactome$Description[3], geneSetID = 3)
# ggsave("gsea_reactome/05_GSEA_Reactome_top1.tiff", plot = reac_gs1, width = 6, height = 5, dpi = 400, compression = "lzw")
# ggsave("gsea_reactome/06_GSEA_Reactome_top2.tiff", plot = reac_gs2, width = 6, height = 5, dpi = 400, compression = "lzw")
# ggsave("gsea_reactome/07_GSEA_Reactome_top3.tiff", plot = reac_gs3, width = 6, height = 5, dpi = 400, compression = "lzw")

reac_gsea1 <- gseaplot2(gsea_reactome, title = gsea_reactome$Description[1], geneSetID = 1)
reac_gsea2 <- gseaplot2(gsea_reactome, title = gsea_reactome$Description[2], geneSetID = 2)
reac_gsea3 <- gseaplot2(gsea_reactome, title = gsea_reactome$Description[3], geneSetID = 3)
# ggsave("gsea_reactome/08_GSEA_Reactome_top1_v2.tiff", plot = reac_gsea1, width = 6, height = 5, dpi = 400, compression = "lzw")
# ggsave("gsea_reactome/09_GSEA_Reactome_top2_v2.tiff", plot = reac_gsea2, width = 6, height = 5, dpi = 400, compression = "lzw")
# ggsave("gsea_reactome/10_GSEA_Reactome_top3_v2.tiff", plot = reac_gsea3, width = 6, height = 5, dpi = 400, compression = "lzw")

reac_gsea5 <- gseaplot2(gsea_reactome, geneSetID = 1:5,
                        color = c("yellow3", "darkblue", "darkolivegreen", "chocolate", "seagreen"))
# ggsave("gsea_reactome/11_GSEA_Reactome_top5_rank.tiff", plot = reac_gsea5, width = 7, height = 5, dpi = 400, compression = "lzw")

# --- 4.5: Network plot ---
reac_cnet <- cnetplot(gsea_reactome, showCategory = 5, categorySize = "pvalue",
                      color.params = list(foldChange = gene_list))
# ggsave("gsea_reactome/12_GSEA_Reactome_cnetplot.tiff", plot = reac_cnet, width = 12, height = 10, dpi = 400, compression = "lzw")

# --- 4.6: Enrichment map ---
gsea_reac_sim <- pairwise_termsim(gsea_reactome)
reac_emap     <- emapplot(gsea_reac_sim, showCategory = 10)
# ggsave("gsea_reactome/13_GSEA_Reactome_emapplot.tiff", plot = reac_emap, width = 10, height = 7, dpi = 400, compression = "lzw")
