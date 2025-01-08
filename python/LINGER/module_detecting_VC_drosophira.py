import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_1samp
import sys
sys.path.append('/data2/duren_lab/cham/cocain/scripts/python/LINGER/')

#from driver_function import *  # Import specific items

Datadir='/data2/duren_lab/cham/cocain/LINGER/'
outdir=Datadir + "output/"#output dir
metadata = pd.read_csv((Datadir + "meta_data2.csv"), index_col=0)
TG_pseudobulk = pd.read_csv((Datadir + "data/TG_pseudobulk.tsv"), index_col=0)
#Gene_score = pd.read_csv("/project/zduren/durenlab/palmetto/cham/Heroin/ori_data/GENE_score_uniq.tsv" ,sep = '\t', index_col=0)
K_list=[7,8,9]
p_list=[0.6,0.7]

trans_reg=pd.read_csv(outdir+'cell_population_trans_regulatory.txt',sep='\t',index_col=0)

DEG_xlsx = Datadir + "combined_anovas.xlsx"
Gene_score_dict = pd.read_excel(DEG_xlsx, sheet_name=None)
score = 'logFC'
score = 'adj.P.Val.Between'
Gene_score_df = make_Gene_score_df(DEG_dir,score)


save_dir = Datadir + 'Module_result_cham/'
score = 'Pr..F.'
compaire = "Line:Treatment\."
compaire = "Line:Treatment:Sex\."

Module_trans_cham_drosophira(trans_reg,metadata,Gene_score_dict,compaire,score,K_list,p_list,save_dir,GWASfile=None,GWAS_score=None,simple=True)
