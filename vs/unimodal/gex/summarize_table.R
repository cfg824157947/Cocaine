library(dplyr)
library(readr)

process_file <- function(file) {
  # Extract cluster number from filename
  cluster_num <- as.numeric(gsub(".*cluster_(\\d+).*", "\\1", file))
  
  # Read raw lines to fix header
  raw_lines <- readLines(file)
  header <- strsplit(raw_lines[1], "\t")[[1]]
  data_line <- strsplit(raw_lines[2], "\t")[[1]]

  # Check if the first value in data row is actually the gene ID or numeric
  header <- c("GeneID", header)

  # Read again with corrected header
  df <- read_delim(file, delim = "\t", col_names = header, skip = 1, col_types = cols())

  # Select FDR columns
  fdr_cols <- c("Treatment_FDR", "Line:Treatment_FDR", "Treatment:Sex_FDR", "Line:Treatment:Sex_FDR")
  df_selected <- df %>% select(all_of(fdr_cols))
  
  # Count FDR < 0.05
  summary_counts <- colSums(df_selected < 0.05, na.rm = TRUE)
  
  # Return result
  return(data.frame(Cluster = cluster_num, t(summary_counts)))
}

# Get all relevant files
files <- list.files(pattern = "cluster_\\d+_seurat_avglogCPM.txt")

# Process all files and combine results
summary_table <- bind_rows(lapply(files, process_file))

# Print summary table
print(summary_table)

# Save summary table to a file
write_delim(summary_table, "summary_fdr_results.txt", delim = "\t")

