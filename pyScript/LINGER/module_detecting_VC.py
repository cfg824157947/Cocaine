import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_1samp
import sys
from scipy import stats
import numpy as np


import numpy as np
from scipy import stats
from sklearn.preprocessing import quantile_transform
from sklearn.decomposition import NMF
#sys.path.append('/project/zduren/durenlab/palmetto/cham/Heroin/script/pyScript/LINGER')

#from driver_function import *  # Import specific items

K_list=[7,8,9]
p_list=[0.7,0.8,0.9]


Datadir='/data2/duren_lab/cham/cocain/LINGER/'
outdir=Datadir + "output/"#output dir
driver_result_dir = Datadir + "driver_results/"
metadata = pd.read_csv((Datadir + "meta_data2.csv"), index_col=0)
#metadata['group'] = (metadata['group'] == 'HEROIN').astype(int)
TG_pseudobulk = pd.read_csv((Datadir + "data/TG_pseudobulk.tsv"), index_col=0)
K=5
adjust_method='bonferroni'
corr_method='pearsonr'


GRN='trans_regulatory'
trans_reg=pd.read_csv(outdir+'cell_population_'+GRN+'.txt',sep='\t',index_col=0)

column_sums=TG_pseudobulk.T.sum()
gene_keep=column_sums[column_sums>10].index
trans_reg=trans_reg[trans_reg.index.isin(gene_keep)]
#TG_pseudobulk=TG_pseudobulk_all
TFset=trans_reg.columns
TGset=trans_reg.index
celltype_list = metadata['celltype'].unique()
#TG_pseudobulk=TG_pseudobulk/TG_pseudobulk.mean()
#idx = [s[:3] != "MIR" for s in TG_pseudobulk.index]
#TG_pseudobulk=TG_pseudobulk.loc[TG_pseudobulk.index[idx]]
R1 = stats.zscore(trans_reg,1);R1[np.isnan(R1)] = 0.0
R2 = stats.zscore(trans_reg,0);R2[np.isnan(R2)] = 0.0
#R1[R1<0]=0;R2any[R2<0]=0
#X=TG_pseudobulk.loc[TGset]
#Exp=quantile_transform(np.log2(X + 1), n_quantiles=10, random_state=0, copy=True)
#Exp=quantile_transform(np.log2(TG_pseudobulk + 1), n_quantiles=10, random_state=0, copy=True)
Exp=pd.DataFrame(TG_pseudobulk,index=TG_pseudobulk.index,columns=TG_pseudobulk.columns)
Z=R1+R2
Z[Z<0]=0  




idx = select_group(metadata,1,"High")
df = make_logFC_df(TG_pseudobulk[idx],metadata.loc[idx],celltype_list)

male_list = metadata['male'].unique()
prefer_list = metadata['prefer'].unique()
Number_of_Significant_TF_df = pd.DataFrame(index = celltype_list)
sns.set(font_scale=1.3)


nmf = NMF(n_components=K, init='random', random_state=0, max_iter=1000)
W = nmf.fit_transform(Z)
H = nmf.components_
W_df=pd.DataFrame(W,index=trans_reg.index,columns=['M'+str(i+1) for i in range(K)])



[S_TG,W2]=assignLabel(W,p);
[S_TF,H2]=assignLabel(H.T,p);
#Exp_mean=stats.zscore(Exp_TG.T).T.groupby(S_TG).mean()
S_TG=pd.DataFrame(S_TG,index=TGset,columns=['Module'])
S_TF=pd.DataFrame(S_TF,index=TFset,columns=['Module'])
Exp_TF=Exp.loc[TFset]
Exp_TG=Exp.loc[TGset]






DEG_dir = "/project/zduren/durenlab/palmetto/cham/Heroin/analyze/NInteResult/RNA/PROTEIN/3/Gene_score/"
score = 'logFC'
score = 'adj.P.Val.Between'
Gene_score_df = make_Gene_score_df(DEG_dir,score)


save_dir = '/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/Module_result_cham/'


Module_trans_cham(trans_reg,metadata,Gene_score_df,K_list,p_list,save_dir,GWASfile=None,GWAS_score=None,simple=True)
