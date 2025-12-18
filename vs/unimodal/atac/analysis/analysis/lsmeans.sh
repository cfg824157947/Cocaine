#!/bin/bash
#
#SBATCH --job-name=CS_filt_atac
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute,gen-mk-compute-1,bigmem,zenyatta,gen-kw-compute-1
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --output=/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/logs/CS_filt_norm.%j.out
#SBATCH --error=/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/logs/CS_filt_norm.%j.err
#SBATCH --mail-type=all
#SBATCH --mail-user=vshanka@clemson.edu

# Get the cluster number from the command-line argument
CLUSTER_NUMBER=$1

# Check if a cluster number was provided
if [ -z "$CLUSTER_NUMBER" ]; then
    echo "Usage: sbatch cs_filter_dor.sh <cluster_number>"
    exit 1
fi

ml R/4.1.2

# Run the R script with the cluster number
#mkdir -p /data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/R/seurat_avglognorm
cd /data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/R/seurat_avglognorm
Rscript /data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/ATAC/R/comp_peaks2DOR/scripts/comp_lsmeans_all_peaks.R $CLUSTER_NUMBER