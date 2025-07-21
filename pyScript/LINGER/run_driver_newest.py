import sys
#sys.path.append('/data2/duren_lab/cham/cocain/scripts/python/LINGER')
#from driver_function import *  # Import specific items
import pandas as pd
import numpy as np   
from statsmodels.stats.multitest import multipletests
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import scipy.stats as stats
from sklearn.decomposition import NMF

#TG_pseudobulk_all=pd.read_csv('data/TG_pseudobulk.tsv',index_col=0,header=0)
Datadir='/data2/duren_lab/cham/cocain/LINGER/'
outdir=Datadir + "output/"#output dir
#driver_result_dir = Datadir + "driver_results/"
driver_result_dir = Datadir + "LINE_result_cham/driver_results_V2/"
metadata = pd.read_csv((Datadir + "meta_data2.csv"), index_col=0)
#metadata['group'] = (metadata['group'] == 'HEROIN').astype(int)
TG_pseudobulk = pd.read_csv((Datadir + "data/TG_pseudobulk.tsv"), index_col=0)
K=25
adjust_method='bonferroni'
corr_method='pearsonr'
#RE_pseudobulk_all=pd.read_csv((Datadir + 'data/RE_pseudobulk.tsv'),index_col=0,header=0)
#TG_pseudobulk_all=pd.read_csv((Datadir + 'data/TG_pseudobulk.tsv'),index_col=0,header=0)


GRN='trans_regulatory'
reg=pd.read_csv(outdir+'cell_population_'+GRN+'.txt',sep='\t',index_col=0)
celltype_list = metadata['celltype'].unique()
metadata['sex'] = metadata['male'].astype(str)
metadata.loc[(metadata.sex == '1'),'sex'] = 'Male'
metadata.loc[(metadata.sex == '0'),'sex'] = 'Female'
metadata['for_group'] = metadata['sex'].astype(str) + '_' + metadata['Line']
for_group_list = metadata['for_group'].unique()
Number_of_Significant_TF_df = pd.DataFrame(index = celltype_list)
sns.set(font_scale=1.3)
melt_df_dict = {}
Write = True
for group in for_group_list:
    idx = select_group_by_key(metadata, 'for_group', group)
    #idx = select_group(metadata,male,Line=prefer)
    df = make_logFC_df(TG_pseudobulk[idx],metadata.loc[idx],celltype_list)
    Gene_score = df.loc[:,(df.isna().sum() == 0)]
    C,P,Q = driver_score_cham_new(reg, Gene_score)
    result_dict ={'C':C, 'P':P, 'Q':Q}
    # C.keys() dict_keys(['original', 'mean', 'std', 'zscore'])
    for key in result_dict.keys():
        for sub_key in result_dict[key].keys():
            result_file = driver_result_dir + key + "_" + sub_key + "_" + group + ".csv"
            if Write:
                result_dict[key][sub_key].to_csv(result_file)


    if Write:
        row_linkage = linkage(result_dict['C']['zscore'], method='average', metric='euclidean')  # Row clustering
        col_linkage = linkage(result_dict['C']['zscore'].T, method='average', metric='euclidean')  # Column clustering
        row_clusters = pd.Series(fcluster(row_linkage, t=5, criterion='maxclust'))
        row_colors = row_clusters.map({i: f"C{i}" for i in np.unique(row_clusters)})
        # Plot a heatmap with clustering on both rows and columns
        clustermap = sns.clustermap(
            result_dict['C']['zscore'],
            row_linkage=row_linkage,
            col_linkage=col_linkage,
            cmap='viridis',     # Color map for the heatmap
            figsize=(15, 15),
            cbar_pos=(1.05, 0.2, 0.03, 0.7)  # Extended margin for color bar
            )
        title = "C_zscore_" + group 
        fig_file = driver_result_dir + title + ".png"
        plt.title(title)
        clustermap.savefig(fig_file, dpi=300)
        plt.clf()

        clustermap = sns.clustermap(
            result_dict['Q']['original'],
            row_linkage=row_linkage,
            col_linkage=col_linkage,
            cmap='viridis',     # Color map for the heatmap
            figsize=(15, 15),
            cbar_pos=(1.05, 0.2, 0.03, 0.7)  # Extended margin for color bar
            )
        title = "Q_original_" + group
        plt.title(title)
        fig_file = driver_result_dir + title + ".png"
        clustermap.savefig(fig_file, dpi=300)
        plt.clf()

    df_C = result_dict['C']['original'].reset_index().melt(
        id_vars='index',
        var_name='CellType',
        value_name='C_value'
    )
    df_CZ = result_dict['C']['zscore'].reset_index().melt(
        id_vars='index',
        var_name='CellType',
        value_name='C_zscore'
    )
    #df_Q = (-np.log10(Q_result_RNA)).reset_index().melt(
    df_Q = (result_dict['Q']['original']).reset_index().melt(
        id_vars='index',
        var_name='CellType',
        value_name='Q_value'
    )


    # 2. Merge on gene name (index) and cell type
    merged_melt_df = pd.merge(df_C, df_Q, on=['index','CellType'])
    merged_melt_df = pd.merge(df_CZ, merged_melt_df, on=['index','CellType'])
    result_file = driver_result_dir + "melt_" + group + ".csv"
    merged_melt_df['Cabs'] = merged_melt_df['C_value'].abs()
    merged_melt_df['group'] = group
    melt_df_dict[group] = merged_melt_df
    if Write:
        merged_melt_df.to_csv(result_file)


#    up_down_list = ['up','down']
#    Q_threshold = 0.01
#    C_threshold = 0.1
#    for up_down in up_down_list:
#        if up_down == 'up':
#            filter_idx = (merged_melt_df['Q_value'] < 0.01) & (merged_melt_df['C_value'] > C_threshold)
#            col_name = sex + "_" + prefer + "_up_C_gt_" + str(C_threshold) + "_Q_lt_" + str(Q_threshold)
#        else:
#            filter_idx = (merged_melt_df['Q_value'] < 0.01) & (merged_melt_df['C_value'] < -C_threshold)
#            col_name = sex + "_" + prefer + "_down_C_lt_" + str(C_threshold) + "_Q_lt_" + str(Q_threshold)
#        count_df = merged_melt_df[filter_idx]['CellType'].value_counts()
#        Number_of_Significant_TF_df.loc[count_df.index, col_name] = count_df


melt_df = pd.concat(melt_df_dict)
melt_df.index = melt_df.group + '_' + melt_df['index']
if Write:
    result_file = driver_result_dir + "all_melt_df.csv"
    melt_df.to_csv(result_file)
#    result_file = driver_result_dir + "permu_Number_of_Significant_TF_df_" + score + ".csv"
#    Number_of_Significant_TF_df.fillna(0).astype(int).to_csv(result_file)


