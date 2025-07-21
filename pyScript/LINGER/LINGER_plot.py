import sys
sys.path.append('/data2/duren_lab/cham/cocain/scripts/python/LINGER')

from driver_function import *  # Import specific items


import pandas as pd
import numpy as np   
from statsmodels.stats.multitest import multipletests
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram


driver_score_dir = "/home/cham/Lab/Cocaine/Cocain_document/LINGER/driver_result/"
result_dir =  "/home/cham/Lab/Cocaine/Cocain_document/LINGER/driver_result/"
file_list = [item for item in os.listdir(driver_score_dir) if ('merged_melt_df_') in item]

df_list=[]
df_dict={}

Q_Thre = 0.01
C_Thre = 0.1

for entry in file_list:
    file_path = driver_score_dir + entry
    df1 =pd.read_csv(file_path,index_col=0)
    group = entry.replace('merged_melt_df_','').replace("_Defalut_logFC.csv","")
    df1['group'] = group
    df1['CellType_group'] = df1['CellType']+'_'+df1['group']
    df1['absC'] = df1['C_value'].abs()
    df1['-log10'] = -np.log10(df1['Q_value'])
#    df_list =df_list+ [df1]
    df_dict[group] =df1
    
for group in df_dict:
    plot_file = result_dir + "dot_" + group + ".png"
#    idx = (df_dict[group]['Q_value']<0.05)
    idx = (df_dict[group]['Q_value']<Q_Thre)& (df_dict[group]['absC']>C_Thre)
#    TF_list = df_dict[group][idx]['index'].unique()
#    idx = df_dict[group]['index'].isin(TF_list)
    plt.figure(figsize=(10, 13))
    scatter = sns.scatterplot(
    #        data=df_dot_plot.loc[idx],
#    data=df_dict[group],
    data=df_dict[group].loc[idx],
    y='index',       # gene name
    x='CellType_group',    # cell/tissue name
    #size='-log10',   # color by C_value
    size='absC',   # color by C_value
    hue='C_value',  # size by Q_value
    sizes=(100, 200), # adjust dot-size range if desired
    palette='viridis'
    )
    scatter.set_xticklabels(scatter.get_xticklabels(), rotation=90)
    # Increase space so labels don’t cut off
    plt.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
    plt.title(group)
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300, bbox_inches="tight")
    #plt.savefig("../../dot_plot2.png", dpi=300, bbox_inches="tight")
    plt.clf()



merged_df = pd.concat((df_dict))

number_of_combination = pd.DataFrame()
idx = (df['Q_value']<Q_Thre)& (df['absC']>C_Thre)
number_of_combination['sum'] = df.loc[idx]['group'].value_counts()
idx = (df['Q_value']<Q_Thre)& (df['C_value']>C_Thre)
number_of_combination['C_mt_' + str(C_Thre)] = df.loc[idx]['group'].value_counts()
idx = (df['Q_value']<Q_Thre)& (df['C_value']<-C_Thre)
number_of_combination['C_lt_-'+str(C_Thre)] = df.loc[idx]['group'].value_counts()
number_of_combination = number_of_combination.sort_index()
result_file = result_dir + "number_of_significant_per_group.csv"
number_of_combination.to_csv(result_file) 


# Create TF_list from Module TF

TF_list_file = "/home/cham/Lab/Cocaine/Cocain_document/LINGER/Module_result/K25_p03_TF.csv"
TF_list = pd.read_csv(TF_list_file, index_col = 0)
TF_list['label_name'] = TF_list.index + '_'+ TF_list['Module'].astype(str)
TF_list = TF_list[~(TF_list['Module']==0)]


TF_list = ['Gsc','CG11294','dl']
#df_dot_plot = merged_df[merged_df['index'].isin(TF_list.index)]
df_dot_plot = merged_df[merged_df['index'].isin(TF_list)]
#df_dot_plot['label_name'] = TF_list.loc[df_dot_plot['index']]['label_name'].values
#df_dot_plot['Module'] = TF_list.loc[df_dot_plot['index']]['Module'].values
#df_dot_plot['-log10'] = -np.log10(df_dot_plot['Q_value'])
#df_dot_plot['absC'] = df_dot_plot.C_value.abs()


#group_idx = df_dot_plot['prefer'] == 'High'
idx = (df_dot_plot['Q_value']<Q_Thre)& (df_dot_plot['absC']>C_Thre)
plt.figure(figsize=(20, 10))
scatter = sns.scatterplot(
    #data=df_dot_plot,
    data=df_dot_plot.loc[idx],
    #data=df_dot_plot.loc[idx&group_idx],
    y='index',       # gene name
    #y='label_name',       # gene name
    x='CellType_group',    # cell/tissue name
    size='absC',   # color by C_value
    hue='C_value',  # size by Q_value
    sizes=(200, 300), # adjust dot-size range if desired
    palette='viridis'
)
scatter.set_xticklabels(scatter.get_xticklabels(), rotation=90)
# Increase space so labels don’t cut off
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
plt.tight_layout()
#title = "M3_TF_driver_result"
#title = 'dot_prefer_High'
#title = 'ALL_Module_TF'
title = 'result from Zhana report'
plt.title(title)
plot_file = result_dir + title + ".png"
plt.savefig(plot_file, dpi=300, bbox_inches="tight")
#plt.savefig("../../dot_plot2.png", dpi=300, bbox_inches="tight")
plt.clf()
