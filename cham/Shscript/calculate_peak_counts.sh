#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=300GB
#SBATCH -o RECALL.out
#SBATCH -e RECALL.err
#SBATCH --mail-type=all
#SBATCH --mail-user=fchang@clemson.edu
#SBATCH --time=48:00:00
##ml anaconda3/2022.05-gcc
source /software/spackages/linux-rocky8-x86_64/gcc-9.5.0/anaconda3-2022.05-zyrazrj6uvrtukupqzhaslr63w7hj6in/etc/profile.d/conda.sh
conda activate /home/fchang/miniconda3/envs/macs2_env


echo "RUN cluster"


celltype_fragment_bed_dir="/project/zduren/durenlab/cham/cocain/ATAC/cluster_fragment/"


celltype_peak_dir="/project/zduren/durenlab/palmetto/cham/cocain/Final/Fragment/Fixed_Recall_peaks/cluster_peaks/"
union_peak_file="/project/zduren/durenlab/palmetto/cham/cocain/Final/Fragment/Fixed_Recall_peaks/union_peak_list.bed" 
celltype_peak_counts_dir="/project/zduren/durenlab/palmetto/cham/cocain/Final/Fragment/Fixed_Recall_peaks/celltype_peak_counts/" 
union_peak_counts_file="/project/zduren/durenlab/palmetto/cham/cocain/Final/Fragment/Fixed_Recall_peaks/union_peak_counts.bed" 

mkdir -p ${celltype_peak_counts_dir}

cat ${celltype_peak_dir}/*.narrowPeak | sort -k1,1 -k2,2n | bedtools merge -i - > ${union_peak_file}


for file in ${celltype_fragment_bed_dir}/*bed
do
	celltype=`basename $file .bed`
	echo $celltype
	echo $file
	echo ${celltype_peak_counts_dir}${celltype}".tsv"

	bedtools intersect -a  ${union_peak_file} -b ${file} -wa -wb| cut -f 1-3,7,8  | sort |bedtools groupby -g 1,2,3,4 -c 5 -o sum> ${celltype_peak_counts_dir}${celltype}".tsv"

#	bedtools intersect -a ${union_peak_file} -b ${file} -wa -wb -sorted | cut -f 1-3,7,8|sort | bedtools groupby -g 1-4 -c 5 -o sum > ${celltype_peak_counts_dir}${celltype}".tsv"
done


cat ${celltype_peak_counts_dir}/* | sort -k1,1 -k2,2n | sed 's/\t/:/;s/\t/-/' > ${union_peak_counts_file}



