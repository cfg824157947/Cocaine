import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import os

# Example DataFrame



C_result_RE,P_result_RE,Q_result_RE=driver_score(RE_pseudobulk_all,metadata,GRN,outdir,adjust_method,corr_method)
C_result_RE_r,Q_result_RE_r=driver_result(C_result_RE,Q_result_RE,K)
C_result_RE_r.to_csv('C_result_RE_r.txt',sep='\t')
Q_result_RE_r.to_csv('Q_result_RE_r.txt',sep='\t')

driver_result_dir = "/data2/duren_lab/cham/cocain/LINGER/driver_result"

file_list = [item for item in os.listdir(driver_result_dir) if 'csv' in item]

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

    # Plot a heatmap with clustering on both rows and columns

    clustermap = sns.clustermap(
        df,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cmap='viridis',     # Color map for the heatmap
        figsize=(15, 15),
        cbar_pos=(1.05, 0.2, 0.03, 0.6)  # Extended margin for color bar
    )

    #plt.tight_layout()
    # Save the figure
    fig_file = os.path.join(driver_result_dir, ("same_order" + entry.split('.')[0]))
    plt.title(entry.split('.')[0])
    clustermap.savefig(fig_file, dpi=300)
    plt.clf()
