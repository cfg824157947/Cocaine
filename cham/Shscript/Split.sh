#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input1> <input2> <output_folder>"
    exit 1
fi

# Variables for readability
input1="$1"
input2="$2"
output_folder="$3"

# Ensure input files exist
if [ ! -f "$input1" ] || [ ! -f "$input2" ]; then
    echo "Error: One or both input files do not exist."
    exit 1
fi

# Process files with awk
awk -v output_folder="$output_folder" '
    BEGIN {
        FS=OFS="\t"
    }
    FNR == 1 && NR > 1 {  # Indicates transition from first file to second
        next_file = 1
        print $0
    }
    !next_file {  # Process the first file
        key_value[$1] = $2
    }
    next_file {  # Process the second file
        if (length(key_value[$4]) != 0) {
            output_file = output_folder "/" key_value[$4] ".bed"
            print $0 >> output_file
        }
    }
' "$input1" "$input2"
