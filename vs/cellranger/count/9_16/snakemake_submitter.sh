#!/bin/bash
#
#SBATCH --job-name=SH10X_multiome_09_16
#SBATCH --ntasks=1   
#SBATCH --partition=bigmem
#SBATCH --time=30-00:00:00
#SBATCH --mem=2gb
#SBATCH --output=/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/log/test_output_%j.txt
#SBATCH --error=/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/log/test_error_%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=vshanka@clemson.edu

cd /data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16
#mkdir -p ./{log,logs_slurm}

source /opt/ohpc/pub/Software/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

#--dag | display | dot
#-p -n \

snakemake \
-s SH10X_multiome_09_16_Snakefile \
--profile slurm \
--configfile SH10X_multiome_09_16.yaml \
--latency-wait 120
