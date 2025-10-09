#! bash
set -euo pipefail

# fragment_dir="/project/zduren/durenlab/palmetto/cham/Heroin/Final/Fragment/Recall_Peak_2ver/sample_celltype/sample_celltype_fragment/"
# cut_dir="/project/zduren/durenlab/palmetto/cham/Heroin/Final/Fragment/Recall_Peak_2ver/sample_celltype/sample_celltype_cut/"
cut_dir="/project/zduren/durenlab/palmetto/cham/cocain/Final/Fragment/cut_dir/"
mkdir -p "$cut_dir"






#TSS_dir="/project/zduren/durenlab/palmetto/cham/util_data/TSS/"
TSS_dir="/project/zduren/durenlab/palmetto/cham/cocain/Final/TSS/"

# for_TSS_dir="/project/zduren/durenlab/palmetto/cham/Heroin/Final/Fragment/Recall_Peak_2ver/sample_celltype/TSS_enrichment/"
for_TSS_dir="/project/zduren/durenlab/palmetto/cham/cocain/Final/Fragment/Fixed_Recall_peaks/sample_cluster/TSS_enrichment/"

TSS_center_bed="${TSS_dir}/TSS_dm6_center.sorted.bed"
TSS_flanks_bed="${TSS_dir}/TSS_dm6_flanks.sorted.bed"



mkdir -p "$for_TSS_dir"

for file in ${cut_dir}*.cuts.bed.gz;
do
    celltype=`basename $file .cuts.bed.gz`
    echo $celltype
    tabix -R ${TSS_center_bed} $file \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -S 2G --parallel=8 \
        | bedtools intersect -sorted -a ${TSS_center_bed} -b - -c \
               > ${for_TSS_dir}${celltype}_TSS_center_counts.tsv

    tabix -R ${TSS_flanks_bed} $file \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -S 2G --parallel=8 \
        | bedtools intersect -sorted -a ${TSS_flanks_bed} -b - -c \
               > ${for_TSS_dir}${celltype}_TSS_flanks_counts.tsv

done




#cd $for_TSS_dir

#for file in ${cut_dir}*.cuts.bed.gz;
#do
#    celltype=`basename $file .cuts.bed.gz`
#    echo $celltype
#    bedtools intersect -u -a $file -b ${TSS_dir}/TSS_flanks.bed | cut -f 1-3 | sed 's/\t/:/;s/\t/-/' > ${cut_dir}/${celltype}_overlap_TSS_flanks.bed
#    bedtools intersect -u -a $file -b ${TSS_dir}/TSS_center.bed | cut -f 1-3 | sed 's/\t/:/;s/\t/-/' > ${cut_dir}/${celltype}_overlap_TSS_center.bed
#done


# bedtools intersect -u -a $file -b ${TSS_dir}/TSS_center.bed | cut -f1-3 | sed 's/\t/:/;s/\t/-/' > ${celltype}_overlap_TSS_center.bed
#what is -u mean>
# -u means write the original A entry once if any overlaps found in B
# bedtools intersect -u -a ../celltype_peak/Astro_peaks.narrowPeak -b ../../../../../util_data/TSS/TSS_center.bed|cut 1-3 | sed 's/\t/:/;s/\t/-/'



