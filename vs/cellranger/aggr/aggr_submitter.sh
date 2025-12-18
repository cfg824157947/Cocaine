#!/bin/bash
#
#SBATCH --job-name=multiome_aggr
#SBATCH --ntasks=30   
#SBATCH --partition=bigmem
#SBATCH --time=30-00:00:00
#SBATCH --mem=750gb
#SBATCH --output=/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Aggr/output_%j.txt
#SBATCH --error=/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Aggr/error_%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=vshanka@clemson.edu

cd /data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Aggr

source /opt/ohpc/pub/Software/anaconda3/etc/profile.d/conda.sh
conda activate gcc9_libstdc6
ml cellranger-arc/2.0.2

cellranger-arc aggr --id=cocaine_multiome_aggr_no_norm \
--csv=/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Aggr/cellranger-arc_aggr_sample_list.csv \
--normalize=none \
--reference=/data/databases/d_mel/d_mel_6.22_ensembl/indexed/10X/cellranger-arc/Dmel6_22 \
--nosecondary
