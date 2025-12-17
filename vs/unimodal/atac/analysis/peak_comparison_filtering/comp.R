#!/usr/bin/env Rscript

# ======================
# Load libraries
# ======================
library(GenomicRanges)
library(dplyr)
library(doParallel)

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
NTR_gr_list <- lapply(NTR_peaks, parse_peak)
NTR_gr <- do.call(c, NTR_gr_list)

# ======================
# Load cluster-specific peaks
# ======================
cluster_bed <- paste0("/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/R/cluster_specific_peaks/",
                      cluster_num, "_peaks.bed")
cluster_peaks <- read.table(cluster_bed, col.names = c("chrom", "start", "end"), sep = "\t")
cluster_gr <- GRanges(seqnames = cluster_peaks$chrom,
                      ranges = IRanges(start = cluster_peaks$start, end = cluster_peaks$end))

# ======================
# Load DOR peaks
# ======================
dor_file <- paste0("/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/R/seurat_avglognorm/cluster_",
                   cluster_num, "_seurat_logCPM_filter.txt")
dor_df <- read.table(dor_file, header = TRUE, sep = "\t")
dor_peaks <- row.names(dor_df)
dor_gr_list <- lapply(dor_peaks, parse_peak)
dor_gr <- do.call(c, dor_gr_list)

# ======================
# Overlap analysis
# ======================
DSC <- subsetByOverlaps(dor_gr, cluster_gr)  # DOR peaks overlapping cluster-specific peaks
DOR_Cluster_specific <- length(DSC)
cluster_specific_DOR <- length(subsetByOverlaps(cluster_gr, dor_gr))
NTR_in_DSC <- length(subsetByOverlaps(NTR_gr, DSC))

DCS_per <- DOR_Cluster_specific / length(dor_gr)
SCD_per <- cluster_specific_DOR / length(cluster_gr)

# ======================
# Print summary
# ======================
cat("\n=== Cluster", cluster_num, "Summary ===\n")
cat("DORs: ", length(dor_gr), "\n")
cat("Cluster-specific Peaks: ", length(cluster_gr), "\n")
cat("DOR overlapping with cluster-specific peaks: ", DOR_Cluster_specific, "\n")
cat("% DOR overlap CSP: ", round(DCS_per * 100, 2), "%\n")
cat("Cluster-specific peaks overlapping with DOR: ", cluster_specific_DOR, "\n")
cat("% CSP overlapping DOR: ", round(SCD_per * 100, 2), "%\n")
cat("NTR Peaks in cluster-specific DORs: ", NTR_in_DSC, "\n\n")

# ======================
# Filter original DOR dataframe for overlapping peaks and save
# ======================
# Convert overlapping GRanges to peak identifiers in "chr:start-end" format
overlapping_peaks <- paste0(seqnames(DSC), ":", start(DSC), "-", end(DSC))

# Filter the DOR dataframe
filtered_dor <- dor_df[rownames(dor_df) %in% overlapping_peaks, , drop = FALSE]

# Output file
out_file <- paste0("cluster_", cluster_num, "_cluster_specific_filtered_DOR.txt")
write.table(filtered_dor, file = out_file, sep = "\t", quote = FALSE, col.names = NA)

cat("Filtered DORs written to:", out_file, "\n")