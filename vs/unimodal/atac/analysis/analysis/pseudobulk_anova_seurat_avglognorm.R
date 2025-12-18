#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript pseudobulk_anova.R <cluster_number>")
}

cluster_number <- args[1]

.libPaths("/opt/ohpc/pub/libs/gnu9/R/4.1.2/lib64/R/library")
library(Signac)
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(reshape2)
library(limma)
library(edgeR)

cocaine_seuratObj <- readRDS(file = "/data2/duren_lab/cham/cocain/Final/Fragments/Union_seurat_obj.rds")

DefaultAssay(cocaine_seuratObj) <- "RNA"

Idents(cocaine_seuratObj) <- "cluster"

pseudobulk <- AverageExpression(cocaine_seuratObj, return.seurat = TRUE, assay = "RNA", layer = "counts", group.by = c('ident','sample_name'), verbose = TRUE, normalization.method="LogNormalize")

cluster_pattern <- paste0(cluster_number, "_")
Cluster_data <- as.data.frame(pseudobulk[["RNA"]]$data[, grep(cluster_pattern, colnames(pseudobulk[["RNA"]]$data), value = TRUE)])


Cluster_data <- Cluster_data[rowMeans(Cluster_data) > 0.015625, ]
Cluster_data <- Cluster_data[rowSums(Cluster_data > 0) >= 4, ]

metadata <- strsplit(sub(cluster_pattern, "", colnames(Cluster_data)), "-")
metadata <- as.data.frame(do.call(rbind, metadata), stringsAsFactors = FALSE)
colnames(metadata) <- c("Line", "Treatment", "Sex", "Replicate")
metadata$Line <- as.factor(metadata$Line)
metadata$Treatment <- as.factor(metadata$Treatment)
metadata$Sex <- as.factor(metadata$Sex)

# Construct a sample_id column that matches the one in your covariate file
metadata$sample_id <- paste0(metadata$Line, "_", metadata$Treatment, "_", metadata$Sex, "_", metadata$Replicate)

# ---- NEW SECTION: Add TSS enrichment covariate ----
# Path to your TSS covariate file
tss_file <- "/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/R/cluster_specific_peaks/TSS_enrichment_summary2k.tsv"

# Read covariate file
tss_data <- read.table(tss_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter for the selected cluster
tss_cluster <- tss_data %>% filter(celltype == paste0("cluster", cluster_number))

# Merge TSS values into metadata by sample_id
metadata <- merge(metadata, tss_cluster[, c("sample_id", "TSS_enrichment")], by = "sample_id", all.x = TRUE)

# Rename the covariate for clarity
colnames(metadata)[colnames(metadata) == "TSS_enrichment"] <- "TSS"

# Check for missing covariates
if (any(is.na(metadata$TSS))) {
  warning("Some samples are missing TSS_enrichment values!")
}


run_anova <- function(gene_expression) {
  lm_fit <- lm(gene_expression ~ Line * Treatment * Sex + TSS, data = metadata)
  anova(lm_fit)
}

anova_results <- apply(Cluster_data, 1, run_anova)

variance_results_list <- list()
terms <- c("Line", "Treatment", "Sex", "Line:Treatment", "Line:Sex", "Treatment:Sex", "Line:Treatment:Sex")

for (i in 1:length(anova_results)) {
  anova_table <- anova_results[[i]]
  p_values <- as.numeric(anova_table[terms, "Pr(>F)"])
  variance_results_list[[i]] <- p_values
}

variance_df <- as.data.frame(do.call(rbind, variance_results_list), stringsAsFactors = FALSE)
colnames(variance_df) <- terms
rownames(variance_df) <- rownames(Cluster_data)

variance_df_fdr <- as.data.frame(apply(variance_df, 2, p.adjust, method = "fdr"))
colnames(variance_df_fdr) <- paste(terms, "FDR", sep = "_")

# Frequency of non-zero samples per gene
non_zero_freq <- rowSums(Cluster_data > 0)
final_variance_df <- cbind(variance_df, variance_df_fdr, NonZero_Sample_Count = non_zero_freq)

final_variance_df_FDR <- final_variance_df[final_variance_df$Treatment_FDR < 0.05 | final_variance_df$`Line:Treatment_FDR` < 0.05 | final_variance_df$`Treatment:Sex_FDR` < 0.05 | final_variance_df$`Line:Treatment:Sex_FDR` < 0.05, ]


output_filename <- paste0("cluster_", cluster_number, "_seurat_logCPM.txt")
output_filename_FDR <- paste0("cluster_", cluster_number, "_seurat_logCPM_filter.txt")

write.table(final_variance_df, output_filename, sep="\t", quote=FALSE, row.names=TRUE)
write.table(final_variance_df_FDR, output_filename_FDR, sep="\t", quote=FALSE, row.names=TRUE)

