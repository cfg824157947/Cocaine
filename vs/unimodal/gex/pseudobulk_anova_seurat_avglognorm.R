.libPaths("/opt/ohpc/pub/libs/gnu9/R/4.1.2/lib64/R/library")
library(Signac)
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(reshape2)
library(limma)
library(edgeR)
library(emmeans)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript pseudobulk_anova.R <cluster_number>")
}

cluster_number <- args[1]

cocaine_seuratObj <- readRDS(file = "/data2/duren_lab/cham/cocain/Data/Coc_seuratObj.rds")

DefaultAssay(cocaine_seuratObj) <- "RNA"

Idents(cocaine_seuratObj) <- "cluster"

pseudobulk <- AverageExpression(cocaine_seuratObj, return.seurat = TRUE, assay = "RNA", layer = "counts", group.by = c('ident','sample_name'), verbose = TRUE, normalization.method="LogNormalize")


cluster_pattern <- paste0("^g",cluster_number, "_")
Cluster_data <- as.data.frame(pseudobulk[["RNA"]]$data[, grep(cluster_pattern, colnames(pseudobulk[["RNA"]]$data), value = TRUE)])

#hist(log(rowMeans(Cluster_data),base = 2),breaks=100)
Cluster_data <- Cluster_data[rowMeans(Cluster_data) > 0.001, ]
Cluster_data <- Cluster_data[rowSums(Cluster_data > 0) >= 4, ]
#Cluster_data <- Cluster_data[rowSums(Cluster_data == 0)/dim(Cluster_data)[2] < 0.96, ]

metadata <- strsplit(sub(cluster_pattern, "", colnames(Cluster_data)), "-")
metadata <- as.data.frame(do.call(rbind, metadata), stringsAsFactors = FALSE)
colnames(metadata) <- c("Line", "Treatment", "Sex", "Replicate")
metadata$Line <- as.factor(metadata$Line)
metadata$Treatment <- as.factor(metadata$Treatment)
metadata$Sex <- as.factor(metadata$Sex)


run_anova <- function(gene_expression) {
  lm_fit <- lm(gene_expression ~ Line * Treatment * Sex, data = metadata)
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


output_filename <- paste0("cluster_", cluster_number, "_seurat_avglogCPM.txt")
output_filename_FDR <- paste0("cluster_", cluster_number, "_seurat_avglogCPM_filter.txt")


write.table(final_variance_df, output_filename, sep="\t", quote=FALSE, row.names=TRUE)
write.table(final_variance_df_FDR, output_filename_FDR, sep="\t", quote=FALSE, row.names=TRUE)

sig_genes <- rownames(Cluster_data)

lsmeans_line_treatment_list <- list()
lsmeans_line_treatment_sex_list <- list()

for (gene in sig_genes) {
  gene_expression <- as.numeric(Cluster_data[gene, ])
  lm_fit <- lm(gene_expression ~ Line * Treatment * Sex, data = metadata)
  
  # ---- Treatment × Line ----
  ls_line_treat <- tryCatch(
    as.data.frame(emmeans(lm_fit, ~ Line * Treatment)),
    error = function(e) NULL
  )
  if (!is.null(ls_line_treat)) {
    # back-transform to original scale
    ls_line_treat$emmean <- exp(ls_line_treat$emmean) - 1
    ls_line_treat$Gene <- gene
    lsmeans_line_treatment_list[[gene]] <- ls_line_treat
  }
  
  # ---- Treatment × Line × Sex ----
  ls_line_treat_sex <- tryCatch(
    as.data.frame(emmeans(lm_fit, ~ Line * Treatment * Sex)),
    error = function(e) NULL
  )
  if (!is.null(ls_line_treat_sex)) {
    # back-transform to original scale
    ls_line_treat_sex$emmean <- exp(ls_line_treat_sex$emmean) - 1
    ls_line_treat_sex$Gene <- gene
    lsmeans_line_treatment_sex_list[[gene]] <- ls_line_treat_sex
  }
}

# ---- Treatment × Line ----
if (length(lsmeans_line_treatment_list) > 0) {
  ls_line_treatment_df <- bind_rows(lsmeans_line_treatment_list)
  
  ls_line_treatment_wide <- ls_line_treatment_df %>%
    mutate(Group = paste(Line, Treatment, sep = "_")) %>%
    select(Gene, Group, emmean) %>%
    pivot_wider(names_from = Group, values_from = emmean)
  
  output_lsmeans_LT <- paste0("cluster_", cluster_number, "_lsmeans_Treatment_Line.txt")
  write.table(ls_line_treatment_wide, output_lsmeans_LT, sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---- Treatment × Line × Sex ----
if (length(lsmeans_line_treatment_sex_list) > 0) {
  ls_line_treatment_sex_df <- bind_rows(lsmeans_line_treatment_sex_list)
  
  ls_line_treatment_sex_wide <- ls_line_treatment_sex_df %>%
    mutate(Group = paste(Line, Treatment, Sex, sep = "_")) %>%
    select(Gene, Group, emmean) %>%
    pivot_wider(names_from = Group, values_from = emmean)
  
  output_lsmeans_LTS <- paste0("cluster_", cluster_number, "_lsmeans_Treatment_Line_Sex.txt")
  write.table(ls_line_treatment_sex_wide, output_lsmeans_LTS, sep = "\t", quote = FALSE, row.names = FALSE)
}
