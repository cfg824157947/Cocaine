#!/bin/bash
#
#SBATCH --job-name=pseudobulk_gex
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute,gen-mk-compute-1,bigmem
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/Gene_expression/logs/seurat_norm.%j.out
#SBATCH --error=/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/Gene_expression/logs/seurat_norm.%j.err
#SBATCH --mail-type=all
#SBATCH --mail-user=vshanka@clemson.edu

# Get the cluster number from the command-line argument
CLUSTER_NUMBER=$1

# Check if a cluster number was provided
if [ -z "$CLUSTER_NUMBER" ]; then
    echo "Usage: sbatch run_pseudobulk_anova.sh <cluster_number>"
    exit 1
fi

ml R/4.1.2

# Run the R script with the cluster number
mkdir -p /data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/Gene_expression/R/seurat_avglognorm
cd /data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/Gene_expression/R/seurat_avglognorm
Rscript /data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/Gene_expression/R/pseudobulk_anova_seurat_avglognorm.R $CLUSTER_NUMBER