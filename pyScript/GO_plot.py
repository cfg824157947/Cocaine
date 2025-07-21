import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

GO_dir = "/home/cham/Lab/Cocaine/Module/"
file_list = [item for item in os.listdir(GO_dir) if ('.csv') in item]
for entry in file_list:
    
    GO_file = os.path.join(GO_dir, entry)
    GO_df = pd.read_csv(GO_file)
    GO_df_sorted = GO_df[(GO_df["P value"]<0.01)]
# Sort GO terms by smallest p-value and take top 5 (adjust as needed)


    # Create the GO enrichment dot plot
    plt.figure(figsize=(20, 30))
    scatter = sns.scatterplot(
        x="P value calc log_10",
        y="Gene Set Name",
        size="Count Overlap Gene",
        hue="Fold Enrichment",
        palette="viridis",
        edgecolor="black",
        sizes=(50, 300),  # Adjust dot sizes
        data=GO_df_sorted
    )

    # Formatting
    plt.xlabel("-log10(P value)", fontsize=30)
    plt.ylabel("GO Term", fontsize=30)
    # Increase font size for x and y ticks
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
#    plt.xlabel("-log10(P value)")
#    plt.ylabel("GO Term")
    title = entry.replace('.csv','GO').replace("Cocain_Data_-_","")
#    title = entry.replace('.csv','GO Analysis')
#    plt.title(title)
    plt.title(title, fontsize=30)
    plt.grid(True, linestyle="--", alpha=0.6)

    # Move the legend further to the right
    legend = plt.legend(loc="lower right")
    plt.setp(legend.get_texts(), fontsize=20)  # Adjust legend font size

    # Save and display
    plt.tight_layout()
    result_name = entry.replace('.csv','_GO_result.png')
    result_file = os.path.join(GO_dir,result_name)
    plt.savefig(result_file, dpi=300)
    plt.clf()
