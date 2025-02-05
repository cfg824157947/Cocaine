import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import os

# Example DataFrame



Datadir='/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/'
driver_result_dir = Datadir + "driver_results/"

score = 'logFC'
file_list = [item for item in os.listdir(driver_result_dir) if (score + '.csv') in item]

for entry in file_list:
    driver_result_file = os.path.join(driver_result_dir, entry)
    print(driver_result_file)

#driver_result_file = "/data2/duren_lab/cham/cocain/LINGER/driver_result/C_TF_TG_driver.csv"
    df = pd.read_csv(driver_result_file, index_col=0)
    df = df.fillna(0)

    if 'C_' in entry:
# Perform hierarchical clustering on rows (TFs) and columns (Cell Types)
        row_linkage = linkage(df, method='average', metric='euclidean')  # Row clustering
        col_linkage = linkage(df.T, method='average', metric='euclidean')  # Column clustering
        row_clusters = pd.Series(fcluster(row_linkage, t=5, criterion='maxclust'))
        row_colors = row_clusters.map({i: f"C{i}" for i in np.unique(row_clusters)})
    # Plot a heatmap with clustering on both rows and columns

    clustermap = sns.clustermap(
        df,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cmap='viridis',     # Color map for the heatmap
        figsize=(15, 15),
        cbar_pos=(1.05, 0.2, 0.03, 0.6)  # Extended margin for color bar
    )
    row_order = clustermap.dendrogram_row.reordered_ind
    col_order = clustermap.dendrogram_col.reordered_ind
    df_reordered = df.iloc[row_order, col_order]
    df_file = os.path.join(driver_result_dir, ("reordered" + entry.split('.')[0] + '.csv'))
    df_reordered.to_csv(df_file)
    #plt.tight_layout()
    # Save the figure
    #fig_file = os.path.join(driver_result_dir, ("same_order" + entry.split('.')[0]))
    fig_file = os.path.join(driver_result_dir, ("same_order" + entry.split('.')[0]))
    plt.title(entry.split('.')[0])
    clustermap.savefig(fig_file, dpi=300)
    plt.clf()

