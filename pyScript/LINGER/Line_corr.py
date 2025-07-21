import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df_dict = {}

def df_row_corr(df1, df2):
    overlap = df1.columns.intersection(df2.columns)
    df1 = df1[overlap]
    df2 = df2[overlap]
    df1_vals = df1.values  # shape (n_rows, n_cols)
    df2_vals = df2.values  # shape (n_rows, n_cols)

    row_correlations = []
    for i in range(df1_vals.shape[0]):
        # np.corrcoef returns a 2x2 matrix when given two 1D arrays
        corr_mat = np.corrcoef(df1_vals[i], df2_vals[i])
        row_correlations.append(corr_mat[0, 1])  # off-diagonal element is the correlation

    row_correlations = np.array(row_correlations)
    row_correlations = pd.DataFrame(row_correlations,index=df1.index)
    return row_correlations

def df_dict_row_corr(df_dict,Module_list=None):
    row_correlations = {}
    if Module_list is None:
        Module_list = df_dict[list(df_dict.keys())[0]].index
    Module_dict = dict(zip(Module_list,range(len(Module_list))))
    for key1 in Module_dict.keys():
        Module_dict[key1] = pd.DataFrame(index = df_dict.keys(),columns=df_dict.keys())

    for key1, df1 in df_dict.items():
        row_correlations[key1] = {}
        for key2, df2 in df_dict.items():
            if key1 == key2:
                for module_name in Module_dict.keys():
                    Module_dict[module_name].loc[key1,key2] = 0
                continue
            row_correlations[key1][key2] = df_row_corr(df1, df2)
            for module_name in Module_dict.keys():
                Module_dict[module_name].loc[key1, key2] = row_correlations[key1][key2].loc[module_name].values
    return row_correlations, Module_dict

df_dict ={}
driver_result_dir = './permu/'
file_list = [item for item in os.listdir(driver_result_dir) if ('Module_C') in item and 'csv' in item and 'mean' in item]
#driver_result_dir = './'
#file_list = [item for item in os.listdir(driver_result_dir) if ('Module_C') in item and 'csv' in item]

for entry in file_list:
    # entry = 'Module_C_result_RNA_Male_RP0231_Pr..F..csv'
    file_path = driver_result_dir + entry
    df = pd.read_csv(file_path,index_col=0)
    index = entry.split('_')[4] + '_' + entry.split('_')[5]
    df_dict[index] = df

row_correlations, Module_dict = df_dict_row_corr(df_dict)

for Module_name in Module_dict.keys():
    Module_dict[Module_name] = Module_dict[Module_name].astype(float)
    row_linkage = sns.clustermap(Module_dict[Module_name], method='average', metric='euclidean')
    col_linkage = sns.clustermap(Module_dict[Module_name].T, method='average', metric='euclidean')
    #plt.figure(figsize=(10, 10))
    clustermap = sns.clustermap(
        Module_dict[Module_name],
        row_linkage=row_linkage.dendrogram_row.linkage,
        col_linkage=col_linkage.dendrogram_col.linkage,
        cmap='viridis',  # Color map for the heatmap
       figsize=(15, 15),
        xticklabels=True,
        yticklabels=True,
    )
    # Adjust font size manually
    clustermap.ax_heatmap.set_xticklabels(clustermap.ax_heatmap.get_xticklabels(), fontsize=15, rotation=45)
    clustermap.ax_heatmap.set_yticklabels(clustermap.ax_heatmap.get_yticklabels(), fontsize=15, rotation=0)
    #sns.heatmap(Module_dict[Module_name].astype(float))
    #plt.xticks(fontsize=8,rotation=45)
    #plt.yticks(fontsize=8)
    title = Module_name
    plt.title(title)
    plt.tight_layout()
    result_dir = "/home/cham/Lab/Cocaine/Cocain_document/LINGER/LINE_result_cham/LINE_corr/"
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    result_file = result_dir + Module_name + '.png'
    plt.savefig(result_file)
    
    plt.close()
    print(Module_name + ' is saved')

'''
# correlation based on the TF score
'''

driver_score_file = driver_result_dir + "all_driver_score.csv"
TF_list = pd.read_csv(driver_score_file,index_col=0)
TF_list = TF_list.loc[((TF_list['C_value']>0.15)&(TF_list['Q_value']<0.01))]['index'].unique()



df_dict ={}
driver_result_dir = './'
#file_list = [item for item in os.listdir(driver_result_dir) if ('C_') in item and 'csv' in item and 'Module' not in item]
file_list = [item for item in os.listdir(driver_result_dir) if ('C_original') in item and 'csv' in item and 'Module' not in item]

for entry in file_list:
    # entry = 'Module_C_result_RNA_Male_RP0231_Pr..F..csv'
    file_path = driver_result_dir + entry
    df = pd.read_csv(file_path,index_col=0)
    df = df.loc[TF_list]
    index = entry.split('_')[2] + '_' + entry.split('_')[3]
    index = index.split('.')[0]
    df_dict[index] = df

row_correlations, Module_dict = df_dict_row_corr(df_dict)

for Module_name in Module_dict.keys():
    Module_dict[Module_name] = Module_dict[Module_name].astype(float)
    row_linkage = sns.clustermap(Module_dict[Module_name], method='average', metric='euclidean')
    col_linkage = sns.clustermap(Module_dict[Module_name].T, method='average', metric='euclidean')
    #plt.figure(figsize=(10, 10))
    clustermap = sns.clustermap(
        Module_dict[Module_name],
        row_linkage=row_linkage.dendrogram_row.linkage,
        col_linkage=col_linkage.dendrogram_col.linkage,
        cmap='viridis',  # Color map for the heatmap
       figsize=(15, 15),
        xticklabels=True,
        yticklabels=True,
    )
    # Adjust font size manually
    clustermap.ax_heatmap.set_xticklabels(clustermap.ax_heatmap.get_xticklabels(), fontsize=15, rotation=45)
    clustermap.ax_heatmap.set_yticklabels(clustermap.ax_heatmap.get_yticklabels(), fontsize=15, rotation=0)
    #sns.heatmap(Module_dict[Module_name].astype(float))
    #plt.xticks(fontsize=8,rotation=45)
    #plt.yticks(fontsize=8)
    title = Module_name
    plt.title(title)
    plt.tight_layout()
    result_dir = "/home/cham/Lab/Cocaine/Cocain_document/LINGER/LINE_result_cham/LINE_corr/V2/"
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    result_file = result_dir + Module_name + '.png'
    plt.savefig(result_file)
    plt.clf()
    plt.close()
    print(Module_name + ' is saved')

# TF patern in each LINE

for entry in file_list:
    # entry = 'Module_C_result_RNA_Male_RP0231_Pr..F..csv'
    file_path = driver_result_dir + entry
    df = pd.read_csv(file_path,index_col=0)
    index = entry.split('_')[2] + '_' + entry.split('_')[3]
    index = index.split('.')[0]
    df_dict[index] = df

row_correlations, Module_dict = df_dict_row_corr(df_dict)

MTFs = []
male_keys = [key for key in row_correlations.keys() if "Male" in key]
#for Line_name in row_correlations.keys():
for Line_name in male_keys:
    TF_corr_LINE = pd.concat(row_correlations[Line_name],axis=1)
    MTFs = MTFs + list(TF_corr_LINE[(TF_corr_LINE.abs()>0.5).sum(axis=1).astype(bool)].index)
# Compute correct hierarchical clustering linkage matrices
    row_linkage = linkage(TF_corr_LINE, method="average", metric="euclidean")
    col_linkage = linkage(TF_corr_LINE.T, method="average", metric="euclidean")
#    row_linkage = sns.clustermap(TF_corr_LINE, method='average', metric='euclidean')
#    col_linkage = sns.clustermap(TF_corr_LINE.T, method='average', metric='euclidean')
    #plt.figure(figsize=(10, 10))
    # Generate clustered heatmap
    clustermap = sns.clustermap(
        TF_corr_LINE,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cmap="viridis",
        figsize=(15, 15),
        xticklabels=True,
        yticklabels=False,
    )
    # Adjust font size manually
    clustermap.ax_heatmap.set_xticklabels(clustermap.ax_heatmap.get_xticklabels(), fontsize=15, rotation=45)
    clustermap.ax_heatmap.set_yticklabels(clustermap.ax_heatmap.get_yticklabels(), fontsize=15, rotation=0)
    title = Line_name
    plt.title(title)
    plt.tight_layout()
    result_dir = "/home/cham/Lab/Cocaine/Cocain_document/LINGER/LINE_result_cham/LINE_corr/V2/"
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    result_file = result_dir + "sep_" + Line_name + '.png'
    plt.savefig(result_file)
    plt.clf()
    plt.close()
    print(Line_name + ' is saved')

# correlation based on the permu mean

df_dict ={}
driver_result_dir = './permu/'
file_list = [item for item in os.listdir(driver_result_dir) if ('Module_C') in item and 'csv' in item and 'mean' in item]

for entry in file_list:
    # entry ='Module_C_result_RNA_mean_Male_RP0231_Defalut_logFC.csv'
    file_path = driver_result_dir + entry
    df = pd.read_csv(file_path,index_col=0)
    index = entry.split('_')[5] + '_' + entry.split('_')[6]
    df_dict[index] = df

row_correlations, Module_dict = df_dict_row_corr(df_dict)

for Module_name in Module_dict.keys():
    Module_dict[Module_name] = Module_dict[Module_name].astype(float)
    row_linkage = sns.clustermap(Module_dict[Module_name], method='average', metric='euclidean')
    col_linkage = sns.clustermap(Module_dict[Module_name].T, method='average', metric='euclidean')
    #plt.figure(figsize=(10, 10))
    clustermap = sns.clustermap(
        Module_dict[Module_name],
        row_linkage=row_linkage.dendrogram_row.linkage,
        col_linkage=col_linkage.dendrogram_col.linkage,
        cmap='viridis',  # Color map for the heatmap
       figsize=(15, 15),
        xticklabels=True,
        yticklabels=True,
    )
    # Adjust font size manually
    clustermap.ax_heatmap.set_xticklabels(clustermap.ax_heatmap.get_xticklabels(), fontsize=15, rotation=45)
    clustermap.ax_heatmap.set_yticklabels(clustermap.ax_heatmap.get_yticklabels(), fontsize=15, rotation=0)
    #sns.heatmap(Module_dict[Module_name].astype(float))
    #plt.xticks(fontsize=8,rotation=45)
    #plt.yticks(fontsize=8)
    title = Module_name
    plt.title(title)
    plt.tight_layout()
    result_dir = "/home/cham/Lab/Cocaine/Cocain_document/LINGER/LINE_result_cham/LINE_corr/permu_mean/"
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    result_file = result_dir + Module_name + '.png'
    plt.savefig(result_file)
    
    plt.close()
    print(Module_name + ' is saved')

