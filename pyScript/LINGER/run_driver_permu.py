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
driver_result_dir = Datadir + "LINE_result_cham/driver_results/permu/"
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
idx = select_group(metadata,1,"High")
df = make_logFC_df(TG_pseudobulk[idx],metadata.loc[idx],celltype_list)

male_list = metadata['male'].unique()
prefer_list = metadata['prefer'].unique()
Line_list = metadata['Line'].unique()
Number_of_Significant_TF_df = pd.DataFrame(index = celltype_list)
score = 'Defalut_logFC'
melt_df_dict = {}
C_result_dict_mean= {}
Q_result_dict_mean= {}

C_result_dict_std= {}
Q_result_dict_std= {}
permu_time = 30
sns.set(font_scale=1.3)
for male in male_list:
    if male == 1:
        sex = 'Male'
    else:
        sex = 'Female'
    #for prefer in prefer_list:
    for prefer in Line_list:
        idx = select_group(metadata,male,Line=prefer)
        df = make_logFC_df(TG_pseudobulk[idx],metadata.loc[idx],celltype_list)
        Gene_score = df.loc[:,(df.isna().sum() == 0)]
        C_result_RNA_mean = None
        P_result_RNA_mean = None
        Q_result_RNA_mean = None
        C_result_std = None
        P_result_std = None
        Q_result_std = None
        for iter in range(permu_time):
            permu_reg = net_col_permute(reg)
            C_result_RNA,P_result_RNA,Q_result_RNA=driver_score_cham(permu_reg, Gene_score)
            if C_result_RNA_mean is None:
                C_result_RNA_mean = C_result_RNA
                P_result_RNA_mean = P_result_RNA
                Q_result_RNA_mean = Q_result_RNA
                C_result_std = C_result_RNA * C_result_RNA
                P_result_std = P_result_RNA * P_result_RNA
                Q_result_std = Q_result_RNA * Q_result_RNA
            else:
                C_result_RNA_mean += C_result_RNA
                P_result_RNA_mean += P_result_RNA
                Q_result_RNA_mean += Q_result_RNA
                C_result_std += C_result_RNA * C_result_RNA
                P_result_std += P_result_RNA * P_result_RNA
                Q_result_std += Q_result_RNA * Q_result_RNA

        C_result_RNA_mean = C_result_RNA_mean / permu_time
        C_result_std = C_result_std / permu_time - C_result_RNA_mean * C_result_RNA_mean
        C_result_std = np.sqrt(C_result_std)
        result_file = driver_result_dir + "C_result_RNA_mean_" + sex + "_" + prefer + "_" + score + ".csv"
        C_result_RNA_mean.fillna(0).to_csv(result_file)
        C_result_dict_mean[sex+'_'+prefer] = C_result_RNA_mean
        C_result_dict_std[sex+'_'+prefer] = C_result_std
        result_file = driver_result_dir + "C_result_RNA_std_" + sex + "_" + prefer + "_" + score + ".csv"
        C_result_std.to_csv(result_file)


        P_result_RNA_mean = P_result_RNA_mean / permu_time
        P_result_std = P_result_std / permu_time - P_result_RNA_mean * P_result_RNA_mean
        P_result_std = np.sqrt(P_result_std)
        result_file = driver_result_dir + "P_result_RNA_mean_" + sex + "_" + prefer + "_" + score + ".csv"
        P_result_RNA_mean.fillna(1).to_csv(result_file)
        result_file = driver_result_dir + "P_result_RNA_std_" + sex + "_" + prefer + "_" + score + ".csv"
        P_result_std.to_csv(result_file)


        Q_result_RNA_mean = Q_result_RNA_mean / permu_time
        Q_result_std = Q_result_std / permu_time - Q_result_RNA_mean * Q_result_RNA_mean
        Q_result_std = np.sqrt(Q_result_std)
        result_file = driver_result_dir + "Q_result_RNA_mean_" + sex + "_" + prefer + "_" + score + ".csv"
        Q_result_RNA_mean.fillna(0).to_csv(result_file)
        Q_result_dict_mean[sex+'_'+prefer] = Q_result_RNA_mean
        Q_result_dict_std[sex+'_'+prefer] = Q_result_std
        result_file = driver_result_dir + "Q_result_RNA_std_" + sex + "_" + prefer + "_" + score + ".csv"
        Q_result_std.to_csv(result_file)

        C_result_RNA = C_result_RNA_mean
        Q_result_RNA = Q_result_RNA_mean
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
        result_file = driver_result_dir + "permu_merged_melt_df_" + sex + "_" + prefer + "_" + score + ".csv"
        merged_melt_df.to_csv(result_file)

        merged_melt_df['C_value'] = merged_melt_df['C_value'].abs()

        melt_df_dict[sex + "_" + prefer] = merged_melt_df


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



result_file = driver_result_dir + "permu_Number_of_Significant_TF_df_" + score + ".csv"
Number_of_Significant_TF_df.fillna(0).astype(int).to_csv(result_file)


