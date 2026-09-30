#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bed_dir> <output_folder>"
    exit 1
fi

# Variables for readability
bed_dir="$1"
out_dir="$2"
#output_folder="$3"
for file in `ls $bed_dir| grep bed`
do

	  echo $file
    celltype=`basename $file .bed`
    echo $cluster
	  macs2 callpeak \
	        -t ${bed_dir}${file} \
	        -f BED \
	        -g dm \
	        --nomodel \
	        --keep-dup all \
	        --call-summits \
	        -n ${out_dir}${cluster}

done
