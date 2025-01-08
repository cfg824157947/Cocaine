import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_1samp
import sys
sys.path.append('/project/zduren/durenlab/palmetto/cham/Heroin/script/pyScript/LINGER')

from driver_function import *  # Import specific items

Datadir='/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/'
outdir=Datadir + "output/"#output dir
metadata = pd.read_csv((Datadir + "data/metadata.csv"), index_col=0)
metadata['group'] = (metadata['group'] == 'HEROIN').astype(int)
TG_pseudobulk = pd.read_csv((Datadir + "data/TG_pseudobulk.tsv"), index_col=0)
Gene_score = pd.read_csv("/project/zduren/durenlab/palmetto/cham/Heroin/ori_data/GENE_score_uniq.tsv" ,sep = '\t', index_col=0)
K_list=[7,8,9]
p_list=[0.7,0.8,0.9]

trans_reg=pd.read_csv(outdir+'cell_population_trans_regulatory.txt',sep='\t',index_col=0)

DEG_dir = "/project/zduren/durenlab/palmetto/cham/Heroin/analyze/NInteResult/RNA/PROTEIN/3/Gene_score/"
score = 'logFC'
score = 'adj.P.Val.Between'
Gene_score_df = make_Gene_score_df(DEG_dir,score)


save_dir = '/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/Module_result_cham/'
Module_trans_cham(trans_reg,metadata,Gene_score_df,K_list,p_list,save_dir,GWASfile=None,GWAS_score=None,simple=True)
