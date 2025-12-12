#!/bin/bash

# Generate nodes file for Cytoscape

Input_dir="/home/cham/Lab/Cocaine/Paper_doc/Table/TG_TF_RE_table_Final/"
Output_dir="/home/cham/Lab/Cocaine/Paper_doc/cytoscape/celltype/"



# Generate edges file for Cytoscape
for file in "${Input_dir}"*"_TG_TF_RE_table.csv"
do
    base_name=$(basename "$file" "_TG_TF_RE_table.csv")
    echo "Generating edges file for ${base_name}"
    edge_file="${Output_dir}${base_name}_edges_for_cytoscape.csv"
    echo "Source,Target,Interaction,Weight" > "$edge_file"
    awk -F, 'NR>1{print $3","$2",regulates,"$4}' "$file" >> "$edge_file"

    # Generate node file for Cytoscape
    echo "Generating nodes file for ${base_name}"
    echo "Id,Label,Type" > "${Output_dir}${base_name}_nodes_for_cytoscape.csv"
    awk -F, 'NR>1{print $2","$2",TG"}' "$file" | sort -u >> "${Output_dir}${base_name}_nodes_for_cytoscape.csv"
    awk -F, 'NR>1{print $3","$3",TF"}' "$file" | sort -u >> "${Output_dir}${base_name}_nodes_for_cytoscape.csv"

done

