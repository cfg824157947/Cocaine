#!/usr/bin/env Rscript

# ======================
# Load libraries
# ======================
.libPaths("/opt/ohpc/pub/libs/gnu9/R/4.1.2/lib64/R/library")
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(Seurat)
library(Signac)
library(Matrix)
library(limma)
library(edgeR)
library(emmeans)

# ======================
# Parse command-line argument
# ======================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript comp.R <cluster_number>")
}
cluster_num <- args[1]  # e.g. "44"

# ======================
# Helper function: parse peak strings (chr:start-end)
# ======================
parse_peak <- function(peak) {
  chrom <- sub(":.*", "", peak)
  pos <- sub(".*:", "", peak)
  start <- as.integer(sub("-.*", "", pos))
  end <- as.integer(sub(".*-", "", pos))
  GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end))
}

# ======================
# Load NTR peaks
# ======================
NTR <- "/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/R/comp_peaks2DOR/NTR_R6.txt"
NTR_peaks <- read.table(NTR, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[[1]]
NTR_gr <- do.call(c, lapply(NTR_peaks, parse_peak))

# ======================
# Load cluster-specific peaks
# ======================
cluster_bed <- paste0("/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/R/cluster_specific_peaks/",
                      cluster_num, "_peaks.bed")
cluster_peaks <- read.table(cluster_bed, col.names = c("chrom", "start", "end"), sep = "\t")
cluster_gr <- GRanges(seqnames = cluster_peaks$chrom,
                      ranges = IRanges(start = cluster_peaks$start, end = cluster_peaks$end))

# ======================
# Load Seurat object
# ======================
cocaine_seuratObj <- readRDS("/data2/duren_lab/cham/cocain/Final/Fragments/Union_seurat_obj.rds")
DefaultAssay(cocaine_seuratObj) <- "RNA"
Idents(cocaine_seuratObj) <- "cluster"

# ======================
# Generate pseudobulk expression per cluster and sample
# ======================
pseudobulk <- AverageExpression(
  cocaine_seuratObj, 
  return.seurat = TRUE, 
  assay = "RNA", 
  layer = "counts", 
  group.by = c("ident","sample_name"), 
  verbose = TRUE, 
  normalization.method = "LogNormalize"
)

# ======================
# Extract cluster-specific pseudobulk data
# ======================
cluster_pattern <- paste0("^g", cluster_num, "_")
Cluster_data <- as.data.frame(pseudobulk[["RNA"]]$data[, grep(cluster_pattern, colnames(pseudobulk[["RNA"]]$data), value = TRUE)])

# Filter lowly expressed peaks
Cluster_data <- Cluster_data[rowMeans(Cluster_data) > 0.015625, ]
Cluster_data <- Cluster_data[rowSums(Cluster_data > 0) >= 4, ]

# ======================
# Build metadata
# ======================
metadata <- strsplit(sub(cluster_pattern, "", colnames(Cluster_data)), "-")
metadata <- as.data.frame(do.call(rbind, metadata), stringsAsFactors = FALSE)
colnames(metadata) <- c("Line", "Treatment", "Sex", "Replicate")
metadata$Line <- as.factor(metadata$Line)
metadata$Treatment <- as.factor(metadata$Treatment)
metadata$Sex <- as.factor(metadata$Sex)
metadata$sample_id <- paste0(metadata$Line, "_", metadata$Treatment, "_", metadata$Sex, "_", metadata$Replicate)

# ======================
# Add TSS enrichment covariate
# ======================
tss_file <- "/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/R/cluster_specific_peaks/TSS_enrichment_summary2k.tsv"
tss_data <- read.table(tss_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tss_cluster <- tss_data %>% filter(celltype == paste0("cluster", cluster_num))
metadata <- merge(metadata, tss_cluster[, c("sample_id", "TSS_enrichment")], by = "sample_id", all.x = TRUE)
colnames(metadata)[colnames(metadata) == "TSS_enrichment"] <- "TSS"

if (any(is.na(metadata$TSS))) {
  warning("Some samples are missing TSS_enrichment values!")
}

# ======================
# Select peaks overlapping cluster-specific peaks
# ======================
all_peaks <- rownames(Cluster_data)
all_peaks_gr <- do.call(c, lapply(all_peaks, parse_peak))
overlapping_peaks_gr <- subsetByOverlaps(all_peaks_gr, cluster_gr)
overlapping_peaks <- paste0(seqnames(overlapping_peaks_gr), ":", start(overlapping_peaks_gr), "-", end(overlapping_peaks_gr))

# Only keep peaks present in Cluster_data
sig_genes <- intersect(rownames(Cluster_data), overlapping_peaks)

# ======================
# Summary statistics
# ======================
cat("\n=== Cluster", cluster_num, "Summary ===\n")
cat("Total peaks in Cluster_data:", nrow(Cluster_data), "\n")
cat("Cluster-specific peaks:", length(cluster_gr), "\n")
cat("Peaks overlapping cluster-specific peaks:", length(sig_genes), "\n")
cat("NTR peaks overlapping cluster-specific peaks:", length(subsetByOverlaps(NTR_gr, overlapping_peaks_gr)), "\n\n")

# ======================
# Calculate lsmeans
# ======================
lsmeans_line_treatment_list <- list()
lsmeans_line_treatment_sex_list <- list()

for (gene in sig_genes) {
  gene_expression <- as.numeric(Cluster_data[gene, ])
  lm_fit <- lm(gene_expression ~ Line * Treatment * Sex, data = metadata)
  
  # Treatment × Line
  ls_line_treat <- tryCatch(
    as.data.frame(emmeans(lm_fit, ~ Line * Treatment)),
    error = function(e) NULL
  )
  if (!is.null(ls_line_treat)) {
    ls_line_treat$emmean <- exp(ls_line_treat$emmean) - 1
    ls_line_treat$Gene <- gene
    lsmeans_line_treatment_list[[gene]] <- ls_line_treat
  }
  
  # Treatment × Line × Sex
  ls_line_treat_sex <- tryCatch(
    as.data.frame(emmeans(lm_fit, ~ Line * Treatment * Sex)),
    error = function(e) NULL
  )
  if (!is.null(ls_line_treat_sex)) {
    ls_line_treat_sex$emmean <- exp(ls_line_treat_sex$emmean) - 1
    ls_line_treat_sex$Gene <- gene
    lsmeans_line_treatment_sex_list[[gene]] <- ls_line_treat_sex
  }
}

# ======================
# Save lsmeans results
# ======================
if (length(lsmeans_line_treatment_list) > 0) {
  ls_line_treatment_df <- bind_rows(lsmeans_line_treatment_list)
  ls_line_treatment_wide <- ls_line_treatment_df %>%
    mutate(Group = paste(Line, Treatment, sep = "_")) %>%
    select(Gene, Group, emmean) %>%
    pivot_wider(names_from = Group, values_from = emmean)
  
  output_lsmeans_LT <- paste0("cluster_", cluster_num, "_lsmeans_Treatment_Line.txt")
  write.table(ls_line_treatment_wide, output_lsmeans_LT, sep = "\t", quote = FALSE, row.names = FALSE)
}

if (length(lsmeans_line_treatment_sex_list) > 0) {
  ls_line_treatment_sex_df <- bind_rows(lsmeans_line_treatment_sex_list)
  ls_line_treatment_sex_wide <- ls_line_treatment_sex_df %>%
    mutate(Group = paste(Line, Treatment, Sex, sep = "_")) %>%
    select(Gene, Group, emmean) %>%
    pivot_wider(names_from = Group, values_from = emmean)
  
  output_lsmeans_LTS <- paste0("cluster_", cluster_num, "_lsmeans_Treatment_Line_Sex.txt")
  write.table(ls_line_treatment_sex_wide, output_lsmeans_LTS, sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("lsmeans calculation completed for cluster", cluster_num, "\n")