#!/usr/bin/env Rscript


# Use commandArgs() to get the input directory from the command line
args <- commandArgs(trailingOnly = TRUE)

# Set input directory location, default to current directory if not provided
input_dir <- if (length(args) > 0) args[1] else "."

# FDR columns to extract
fdr_cols <- c("Treatment_FDR", "Line:Treatment_FDR", "Treatment:Sex_FDR", "Line:Treatment:Sex_FDR")

# Store merged results for each FDR column
fdr_results <- list()
for (fdr in fdr_cols) fdr_results[[fdr]] <- NULL

# Loop through cluster files
for (i in 0:47) {
  file_path <- file.path(input_dir, paste0("cluster_", i, "_seurat_logCPM_filter.txt"))

  if (!file.exists(file_path)) {
    warning(paste("Missing:", file_path))
    next
  }

  # Skip if file only has the header line
  if (length(readLines(file_path, warn = FALSE)) <= 1) {
    next
  }
  
  #raw <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE)
  #headers <- as.character(unlist(raw[1, -1]))
  df <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE, skip=1)
  colnames(df) <- c("Gene", "Line" , "Treatment", "Sex", "Line:Treatment",  "Line:Sex",  "Treatment:Sex", "Line:Treatment:Sex",  "Line_FDR",  "Treatment_FDR", "Sex_FDR", "Line:Treatment_FDR",  "Line:Sex_FDR",  "Treatment:Sex_FDR", "Line:Treatment:Sex_FDR")
  #df[df == "<NA>"] <- NA
  df <- df[rowSums(is.na(df)) < 14, ]
  rownames(df) <- df$Gene
  df$Gene <- NULL
  #df <- df[rowSums(!is.na(df)) > 0, ]

  for (fdr in fdr_cols) {
    if (fdr %in% colnames(df)) {
      sub_df <- data.frame(Gene = rownames(df), value = as.numeric(df[[fdr]]), stringsAsFactors = FALSE)
      colnames(sub_df)[2] <- paste0("Cluster_", i)
      if (is.null(fdr_results[[fdr]])) {
        fdr_results[[fdr]] <- sub_df
      } else {
        fdr_results[[fdr]] <- merge(fdr_results[[fdr]], sub_df, by = "Gene", all = TRUE)
      }
    }
  }
}

# Write out the merged results
for (fdr in fdr_cols) {
  out_df <- fdr_results[[fdr]]
  rownames(out_df) <- out_df$Gene
  out_df$Gene <- NULL
  out_df <- t(out_df)  # Cluster as row, genes as columns
  write.csv(out_df, paste0(fdr, "_by_cluster.csv"), quote = FALSE, na = "")
}
