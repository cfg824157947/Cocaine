#import sys
#sys.path.append('/project/zduren/durenlab/palmetto/cham/Heroin/script/pyScript/LINGER')

#from driver_function import *  # Import specific items


import pandas as pd
import numpy as np   
from statsmodels.stats.multitest import multipletests
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
#TG_pseudobulk_all=pd.read_csv('data/TG_pseudobulk.tsv',index_col=0,header=0)
Datadir='/data2/duren_lab/cham/cocain/LINGER/'
outdir=Datadir + "output/"#output dir
driver_result_dir = Datadir + "driver_results/"
metadata = pd.read_csv((Datadir + "meta_data2.csv"), index_col=0)
#metadata['group'] = (metadata['group'] == 'HEROIN').astype(int)
TG_pseudobulk = pd.read_csv((Datadir + "data/TG_pseudobulk.tsv"), index_col=0)
K=5
adjust_method='bonferroni'
corr_method='pearsonr'
#RE_pseudobulk_all=pd.read_csv((Datadir + 'data/RE_pseudobulk.tsv'),index_col=0,header=0)
#TG_pseudobulk_all=pd.read_csv((Datadir + 'data/TG_pseudobulk.tsv'),index_col=0,header=0)


GRN='trans_regulatory'
reg=pd.read_csv(outdir+'cell_population_'+GRN+'.txt',sep='\t',index_col=0)

idx = select_group(metadata,1,"High")
df = make_logFC_df(TG_pseudobulk[idx],metadata.loc[idx],celltype_list)
celltype_list = metadata['celltype'].unique()

male_list = metadata['male'].unique()
prefer_list = metadata['prefer'].unique()
Number_of_Significant_TF_df = pd.DataFrame(index = celltype_list)
sns.set(font_scale=1.3)
for male in male_list:
    if male == 1:
        sex = 'Male'
    else:
        sex = 'Female'
    for prefer in prefer_list:
        idx = select_group(metadata,male,prefer)
        df = make_logFC_df(TG_pseudobulk[idx],metadata.loc[idx],celltype_list)
        Gene_score = df.loc[:,(df.isna().sum() == 0)]
        C_result_RNA,P_result_RNA,Q_result_RNA=driver_score_cham(reg, Gene_score)
        result_file = driver_result_dir + "C_result_RNA_" + sex + "_" + prefer + "_" + score + ".csv"
        C_result_RNA.fillna(0).to_csv(result_file)
        result_file = driver_result_dir + "Q_result_RNA_" + sex + "_" + prefer + "_" + score + ".csv"
        Q_result_RNA.fillna(1).to_csv(result_file)

        row_linkage = linkage(C_result_RNA, method='average', metric='euclidean')  # Row clustering
        col_linkage = linkage(C_result_RNA.T, method='average', metric='euclidean')  # Column clustering
        row_clusters = pd.Series(fcluster(row_linkage, t=5, criterion='maxclust'))
        row_colors = row_clusters.map({i: f"C{i}" for i in np.unique(row_clusters)})

        # Plot a heatmap with clustering on both rows and columns
        clustermap = sns.clustermap(
            C_result_RNA,
            row_linkage=row_linkage,
            col_linkage=col_linkage,
            cmap='viridis',     # Color map for the heatmap
            figsize=(15, 15),
            cbar_pos=(1.05, 0.2, 0.03, 0.7)  # Extended margin for color bar
            )
        title = "C_result_RNA_" + sex + "_" + prefer + "_" + score 
        fig_file = driver_result_dir + title + ".png"
        plt.title(title)
        clustermap.savefig(fig_file, dpi=300)
        plt.clf()

        clustermap = sns.clustermap(
            Q_result_RNA,
            row_linkage=row_linkage,
            col_linkage=col_linkage,
            cmap='viridis',     # Color map for the heatmap
            figsize=(15, 15),
            cbar_pos=(1.05, 0.2, 0.03, 0.7)  # Extended margin for color bar
            )
        title = "Q_result_RNA_" + sex + "_" + prefer + "_" + score 
        plt.title(title)
        fig_file = driver_result_dir + title + ".png"
        clustermap.savefig(fig_file, dpi=300)
        plt.clf()

        df_C = C_result_RNA.reset_index().melt(
            id_vars='index',
            var_name='CellType',
            value_name='C_value'
        )
        #df_Q = (-np.log10(Q_result_RNA)).reset_index().melt(
        df_Q = (Q_result_RNA).reset_index().melt(
            id_vars='index',
            var_name='CellType',
            value_name='Q_value'
        )

        # 2. Merge on gene name (index) and cell type
        merged_melt_df = pd.merge(df_C, df_Q, on=['index','CellType'])
        result_file = driver_result_dir + "merged_melt_df_" + sex + "_" + prefer + "_" + score + ".csv"
        merged_melt_df.to_csv(result_file)

        merged_melt_df['C_value'] = merged_melt_df['C_value'].abs()

        up_down_list = ['up','down']
        Q_threshold = 0.01
        C_threshold = 0.1
        for up_down in up_down_list:
            if up_down == 'up':
                filter_idx = (merged_melt_df['Q_value'] < 0.01) & (merged_melt_df['C_value'] > C_threshold)
                col_name = sex + "_" + prefer + "_up_C_gt_" + str(C_threshold) + "_Q_lt_" + str(Q_threshold)
            else:
                filter_idx = (merged_melt_df['Q_value'] < 0.01) & (merged_melt_df['C_value'] < -C_threshold)
                col_name = sex + "_" + prefer + "_down_C_lt_" + str(C_threshold) + "_Q_lt_" + str(Q_threshold)
            count_df = merged_melt_df[filter_idx]['CellType'].value_counts()
            Number_of_Significant_TF_df.loc[count_df.index, col_name] = count_df



result_file = driver_result_dir + "Number_of_Significant_TF_df_" + score + ".csv"
Number_of_Significant_TF_df.fillna(0).astype(int).to_csv(result_file)






        # 1. Get the top-3 columns for each row (as a list)
        top3_TF_per_celltype = C_result_RNA.T.apply(lambda row: row.nlargest(1).index.tolist(), axis=1)
        # 2. Flatten and convert to a set to remove duplicates
        unique_top_TF = set(col for row_list in top3_TF_per_celltype for col in row_list)
        # 3. Convert back to a list if desired
        TF_list = list(unique_top_TF)

        df_dot_plot = merged_melt_df[merged_melt_df['index'].isin(TF_list)]





        plt.figure(figsize=(20, 40))

        scatter = sns.scatterplot(
            data=df_plot,
            y='index',       # gene name
            x='CellType',    # cell/tissue name
            hue='C_value',   # color by C_value
            size='Q_value',  # size by Q_value
            sizes=(20, 200), # adjust dot-size range if desired
            palette='viridis'
        )
        scatter.set_xticklabels(scatter.get_xticklabels(), rotation=90)
        # Increase space so labels don’t cut off
        plt.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
        plt.tight_layout()
        plt.savefig("./dot_plot.png", dpi=300, bbox_inches="tight")
        plt.clf()





GRN='TF_RE_binding'


idx = select_group(metadata,1,"High")
df = make_logFC_df(TG_pseudobulk[idx],metadata.loc[idx],celltype_list)


C_result_RNA,P_result_RNA,Q_result_RNA=driver_score_cham(reg, Gene_score_df)
C_result_RNA,P_result_RNA,Q_result_RNA=driver_score_cham(reg,adjust_method,corr_method, Gene_score_df)

result_file = driver_result_dir + "C_result_RNA_" + score + ".csv"
C_result_RNA.fillna(0).to_csv(result_file)

row_linkage = linkage(C_result_RNA, method='average', metric='euclidean')  # Row clustering
col_linkage = linkage(C_result_RNA.T, method='average', metric='euclidean')  # Column clustering
row_clusters = pd.Series(fcluster(row_linkage, t=5, criterion='maxclust'))
row_colors = row_clusters.map({i: f"C{i}" for i in np.unique(row_clusters)})

# Plot a heatmap with clustering on both rows and columns
clustermap = sns.clustermap(
    C_result_RNA,
    row_linkage=row_linkage,
    col_linkage=col_linkage,
    cmap='viridis',     # Color map for the heatmap
    figsize=(15, 15),
    cbar_pos=(1.05, 0.2, 0.03, 0.6)  # Extended margin for color bar
    )



result_file = driver_result_dir + "Q_result_RNA_" + score + ".csv"
Q_result_RNA.fillna(1).to_csv(result_file)

clustermap = sns.clustermap(
    Q_result_RNA,
    row_linkage=row_linkage,
    col_linkage=col_linkage,
    cmap='viridis',     # Color map for the heatmap
    figsize=(15, 15),
    cbar_pos=(1.05, 0.2, 0.03, 0.6)  # Extended margin for color bar
    )




GRN='TF_RE_binding'
reg=pd.read_csv(outdir+'cell_population_'+GRN+'.txt',sep='\t',index_col=0)

DAR_dir = "/project/zduren/durenlab/palmetto/cham/Heroin/analyze/NInteResult/ATAC/WHOLE/3/Gene_score/"

score = 'logFC'
Gene_score_df = make_Gene_score_df(DAR_dir,score)
C_result_RNA,P_result_RNA,Q_result_RNA=driver_score_cham(reg, Gene_score_df)

result_file = driver_result_dir + "C_result_RE_" + score + ".csv"
C_result_RNA.fillna(0).to_csv(result_file)
result_file = driver_result_dir + "Q_result_RE_" + score + ".csv"
Q_result_RNA.fillna(1).to_csv(result_file)
