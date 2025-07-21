import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_1samp
import sys
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import numpy as np
from sklearn.decomposition import NMF
import scipy.stats as stats
sys.path.append('/data2/duren_lab/cham/cocain/scripts/python/LINGER/')

#from driver_function import *  # Import specific items

Datadir='/data2/duren_lab/cham/cocain/LINGER/'
outdir=Datadir + "output/"#output dir
metadata = pd.read_csv((Datadir + "meta_data2.csv"), index_col=0)
TG_pseudobulk = pd.read_csv((Datadir + "data/TG_pseudobulk.tsv"), index_col=0)
#Gene_score = pd.read_csv("/project/zduren/durenlab/palmetto/cham/Heroin/ori_data/GENE_score_uniq.tsv" ,sep = '\t', index_col=0)
K_list=[7,8,9]
p_list=[0.6,0.7]
K=25

trans_reg=pd.read_csv(outdir+'cell_population_trans_regulatory.txt',sep='\t',index_col=0)

#DEG_xlsx = Datadir + "combined_anovas.xlsx"
#Gene_score_dict = pd.read_excel(DEG_xlsx, sheet_name=None)
#score = 'logFC'
#score = 'adj.P.Val.Between'
#Gene_score_df = make_Gene_score_df(DEG_dir,score)


save_dir = Datadir + 'Module_result_cham/'
save_dir = Datadir + 'Module_result_nosex_cham/'
save_dir = Datadir + 'LINE_result_cham/'
score = 'Pr..F.'
compaire = "Line:Treatment:Sex\."
compaire = "Line:Treatment\."


print('loading GRN......')
column_sums=TG_pseudobulk.T.sum()
gene_keep=column_sums[column_sums>10].index
trans_reg=trans_reg[trans_reg.index.isin(gene_keep)]
#TG_pseudobulk=TG_pseudobulk_all
TFset=trans_reg.columns
TGset=trans_reg.index
celltype = metadata['celltype'].unique()

R1 = stats.zscore(trans_reg,1);R1[np.isnan(R1)] = 0.0
R2 = stats.zscore(trans_reg,0);R2[np.isnan(R2)] = 0.0
Exp=pd.DataFrame(TG_pseudobulk,index=TG_pseudobulk.index,columns=TG_pseudobulk.columns)
Z=R1+R2
Z[Z<0]=0  





# Generate Module
print('identify modules......')
nmf = NMF(n_components=K, init='random', random_state=0, max_iter=1000)
W = nmf.fit_transform(Z)
H = nmf.components_
W_df=pd.DataFrame(W,index=trans_reg.index,columns=['M'+str(i+1) for i in range(K)])


#driver_result_dir = Datadir + "driver_results/"
driver_result_dir = Datadir + "LINE_result_cham/driver_results/permu/"
male_list = metadata['male'].unique()
prefer_list = metadata['prefer'].unique()
Line_list = metadata['Line'].unique()
celltype_list = metadata['celltype'].unique()
K=25
C_result_dict_mean= {}
Q_result_dict_mean= {}

C_result_dict_std= {}
Q_result_dict_std= {}
permu_time = 30
for male in male_list:
    if male == 1:
        sex = 'Male'
    else:
        sex = 'Female'
    #    for prefer in prefer_list:
    for prefer in Line_list:
#        idx = select_group(metadata,male,prefer)
        idx = select_group(metadata,male=male,Line=prefer)
        df = make_logFC_df(TG_pseudobulk[idx],metadata.loc[idx],celltype_list)
        Gene_score = df.loc[:,(df.isna().sum() == 0)]
        C_result_RNA_mean = None
        P_result_RNA_mean = None
        Q_result_RNA_mean = None
        C_result_std = None
        P_result_std = None
        Q_result_std = None
        for iter in range(permu_time):
            permu_W_df = net_col_permute(W_df)
            C_result_RNA,P_result_RNA,Q_result_RNA=driver_score_cham(permu_W_df, Gene_score)
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
        result_file = driver_result_dir + "Module_C_result_RNA_mean_" + sex + "_" + prefer + "_" + score + ".csv"
        C_result_RNA_mean.fillna(0).to_csv(result_file)
        C_result_dict_mean[sex+'_'+prefer] = C_result_RNA_mean
        C_result_dict_std[sex+'_'+prefer] = C_result_std
        result_file = driver_result_dir + "Module_C_result_RNA_std_" + sex + "_" + prefer + "_" + score + ".csv"
        C_result_std.to_csv(result_file)


        P_result_RNA_mean = P_result_RNA_mean / permu_time
        P_result_std = P_result_std / permu_time - P_result_RNA_mean * P_result_RNA_mean
        P_result_std = np.sqrt(P_result_std)
        result_file = driver_result_dir + "Module_P_result_RNA_mean_" + sex + "_" + prefer + "_" + score + ".csv"
        P_result_RNA_mean.fillna(1).to_csv(result_file)
        result_file = driver_result_dir + "Module_P_result_RNA_std_" + sex + "_" + prefer + "_" + score + ".csv"
        P_result_std.to_csv(result_file)



        #result_file = driver_result_dir + "Module_P_result_RNA_mean_" + sex + "_" + prefer + "_" + score + ".csv"
        #P_result_RNA.fillna(1).to_csv(result_file)

        Q_result_RNA_mean = Q_result_RNA_mean / permu_time
        Q_result_std = Q_result_std / permu_time - Q_result_RNA_mean * Q_result_RNA_mean
        Q_result_std = np.sqrt(Q_result_std)
        result_file = driver_result_dir + "Module_Q_result_RNA_mean_" + sex + "_" + prefer + "_" + score + ".csv"
        Q_result_RNA_mean.fillna(0).to_csv(result_file)
        Q_result_dict_mean[sex+'_'+prefer] = Q_result_RNA_mean
        Q_result_dict_std[sex+'_'+prefer] = Q_result_std
        result_file = driver_result_dir + "Module_Q_result_RNA_std_" + sex + "_" + prefer + "_" + score + ".csv"
        Q_result_std.to_csv(result_file)

        C_result_RNA = C_result_RNA_mean
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
        title = "Module_C_result_RNA_mean_" + sex + "_" + prefer + "_" + score 
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
        title = "Module_Q_result_RNA_mean_" + sex + "_" + prefer + "_" + score 
        plt.title(title)
        fig_file = driver_result_dir + title + ".png"
        clustermap.savefig(fig_file, dpi=300)
        plt.clf()


[S_TG,W2]=assignLabel(W,p);
[S_TF,H2]=assignLabel(H.T,p);
#Exp_mean=stats.zscore(Exp_TG.T).T.groupby(S_TG).mean()
S_TG=pd.DataFrame(S_TG,index=TGset,columns=['Module'])
S_TF=pd.DataFrame(S_TF,index=TFset,columns=['Module'])
Exp_TF=Exp.loc[TFset]
Exp_TG=Exp.loc[TGset]




#Module_trans_cham_drosophira(trans_reg,metadata,Gene_score_dict,compaire,score,K_list,p_list,save_dir,GWASfile=None,GWAS_score=None,simple=True)
