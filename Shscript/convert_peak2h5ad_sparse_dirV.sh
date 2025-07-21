#!/bin/bash

# Define base directories
bed_dir="/project/zduren/durenlab/palmetto/cham/cocain/Final/Fragment/Fixed_Recall_peaks/sc_cluster_peak_counts/bed_dir"
h5ad_dir="/project/zduren/durenlab/palmetto/cham/cocain/Final/Fragment/Fixed_Recall_peaks/sc_cluster_peak_counts/h5ad_dir"
sparse_dir="/project/zduren/durenlab/palmetto/cham/cocain/Final/Fragment/Fixed_Recall_peaks/sc_cluster_peak_counts/sparse_matrix_dir"
script_path="/project/zduren/durenlab/palmetto/cham/cocain/Final/git_handle/Cocaine/pyScript/convert_peak_to_adata.py"  # ⬅️ Update this to your actual script location

# Loop through each TSV file
for tsv_file in "$bed_dir"/*.tsv; do
    base_name=$(basename "$tsv_file" .tsv)
    
    tsv_path="$tsv_file"
    h5ad_path="${h5ad_dir}/${base_name}.h5ad"
    outdir="${sparse_dir}/${base_name}_sparse_matrix"

    echo "Processing $base_name ..."
    python "$script_path" \
           --tsv "$tsv_path" \
           --h5ad "$h5ad_path" \
           --outdir "$outdir" \
           --feature_type peak
done

echo "All files processed."
