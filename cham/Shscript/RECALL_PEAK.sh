#!/bin/bash
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 12-00:00
#SBATCH --mem=64GB
#SBATCH -p gen-mk-compute-1
#SBATCH -o RECALL.out
#SBATCH -e RECALL.err
#SBATCH --mail-type=all
#SBATCH --mail-user=fchang@clemson.edu


echo "RUN cluster"

label_file="/data2/duren_lab/cham/cocain/ATAC/condition_label.tsv"
fragments_file="/data2/duren_lab/cham/cocain/ATAC/atac_fragments.tsv"
bed_dir="/data2/duren_lab/cham/cocain/ATAC/cluster_fragment/"
peak_dir="/data2/duren_lab/cham/cocain/ATAC/peak_result/"
peak_count_dir="/data2/duren_lab/cham/cocain/ATAC/peak_result/peak_counts/"
sh /data2/duren_lab/cham/cocain/ATAC/Split.sh $label_file $fragments_file $bed_dir
sh /data2/duren_lab/cham/cocain/ATAC/umacs2.sh $bed_dir $peak_dir
sh /project/zduren/durenlab/palmetto/cham/Heroin/ATAC/new/calculate_peak_counts.sh ${bed_dir} ${peak_dir} ${peak_count_dir}
