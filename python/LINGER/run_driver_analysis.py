import sys
sys.path.append('/project/zduren/durenlab/palmetto/cham/Heroin/script/pyScript/LINGER')

from driver_function import *  # Import specific items


import pandas as pd
import numpy as np   
from statsmodels.stats.multitest import multipletests
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
RE_pseudobulk_all=pd.read_csv('data/RE_pseudobulk.tsv',index_col=0,header=0)
TG_pseudobulk_all=pd.read_csv('data/TG_pseudobulk.tsv',index_col=0,header=0)
Datadir='/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/'
outdir=Datadir + "output/"#output dir
driver_result_dir = Datadir + "driver_results/"
K=5
adjust_method='bonferroni'
corr_method='pearsonr'


DEG_dir = "/project/zduren/durenlab/palmetto/cham/Heroin/analyze/NInteResult/RNA/PROTEIN/3/Gene_score/"
GRN='trans_regulatory'
reg=pd.read_csv(outdir+'cell_population_'+GRN+'.txt',sep='\t',index_col=0)

score = 'logFC'
score = 'adj.P.Val.Between'
Gene_score_df = make_Gene_score_df(DEG_dir,score)

C_result_RNA,P_result_RNA,Q_result_RNA=driver_score_cham(reg,adjust_method,corr_method, Gene_score_df)

result_file = driver_result_dir + "C_result_RNA_" + score + ".csv"
C_result_RNA.fillna(0).to_csv(result_file)
result_file = driver_result_dir + "Q_result_RNA_" + score + ".csv"
Q_result_RNA.fillna(1).to_csv(result_file)


GRN='TF_RE_binding'
reg=pd.read_csv(outdir+'cell_population_'+GRN+'.txt',sep='\t',index_col=0)

DAR_dir = "/project/zduren/durenlab/palmetto/cham/Heroin/analyze/NInteResult/ATAC/WHOLE/3/Gene_score/"

score = 'logFC'
Gene_score_df = make_Gene_score_df(DAR_dir,score)
C_result_RNA,P_result_RNA,Q_result_RNA=driver_score_cham(reg,adjust_method,corr_method, Gene_score_df)

result_file = driver_result_dir + "C_result_RE_" + score + ".csv"
C_result_RNA.fillna(0).to_csv(result_file)
result_file = driver_result_dir + "Q_result_RE_" + score + ".csv"
Q_result_RNA.fillna(1).to_csv(result_file)
