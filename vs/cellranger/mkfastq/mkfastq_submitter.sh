#!/bin/bash
#
#SBATCH --job-name=SH10X_ATAC_1_8
#SBATCH --ntasks=1   
#SBATCH --partition=compute
#SBATCH --time=30-00:00:00
#SBATCH --mem=2gb
#SBATCH --output=/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/FASTQS/SH10X_ATAC_1_8/output_%j.txt
#SBATCH --error=/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/FASTQS/SH10X_ATAC_1_8/error_%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=vshanka@clemson.edu

cd /data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/FASTQS/SH10X_ATAC_1_8

source /opt/ohpc/pub/Software/anaconda3/etc/profile.d/conda.sh
conda activate gcc9_libstdc6
ml cellranger-arc/2.0.2

cellranger-arc mkfastq --id=SH10X_ATAC_1_8 \
--run=/data/Palmetto_sync/Novaseq/BCL/scCocaine_multiome_08_2022/SH10x_ATAC_1_8/Files \
--csv=ATACSeq_multiome_SH_01-08_corrected.csv \
--jobmode=slurm
