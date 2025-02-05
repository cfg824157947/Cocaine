import LingerGRN.LL_net as LL_net
import LingerGRN.LINGER_tr as LINGER_tr
from LingerGRN.preprocess import *
from LingerGRN.pseudo_bulk import *
import os
import scipy.sparse as sp
import pandas as pd
import scanpy as sc

from LingerGRN.preprocess import *

h5_file = "/project/zduren/durenlab/palmetto/zduren/DrachevaLab/HeroinData/merged/aggregation_withoutnormlize/outs/filtered_feature_bc_matrix.h5"
label_file = "/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/label.csv"

adata = sc.read_10x_h5(h5_file ,gex_only=False)
label = pd.read_csv(label_file)


matrix=adata.X.T
adata.var['gene_ids']=adata.var.index
features=pd.DataFrame(adata.var['gene_ids'].values.tolist(),columns=[1])
features[2]=adata.var['feature_types'].values
barcodes=pd.DataFrame(adata.obs_names,columns=[0])
adata_RNA,adata_ATAC=get_adata(matrix,features,barcodes,label)# adata_RNA and adata_ATAC are scRNA and scATAC

sc.pp.filter_cells(adata_RNA, min_genes=200)
sc.pp.filter_genes(adata_RNA, min_cells=3)
sc.pp.filter_cells(adata_ATAC, min_genes=200)
sc.pp.filter_genes(adata_ATAC, min_cells=3)
selected_barcode=list(set(adata_RNA.obs['barcode'].values)&set(adata_ATAC.obs['barcode'].values))
barcode_idx=pd.DataFrame(range(adata_RNA.shape[0]), index=adata_RNA.obs['barcode'].values)
adata_RNA = adata_RNA[barcode_idx.loc[selected_barcode][0]]
barcode_idx=pd.DataFrame(range(adata_ATAC.shape[0]), index=adata_ATAC.obs['barcode'].values)
adata_ATAC = adata_ATAC[barcode_idx.loc[selected_barcode][0]]


samplelist=list(set(adata_ATAC.obs['sample'].values)) # sample is generated from cell barcode 
tempsample=samplelist[0]
TG_pseudobulk=pd.DataFrame([])
RE_pseudobulk=pd.DataFrame([])
singlepseudobulk = (adata_RNA.obs['sample'].unique().shape[0]*adata_RNA.obs['sample'].unique().shape[0]>100)
for tempsample in samplelist:
    adata_RNAtemp=adata_RNA[adata_RNA.obs['sample']==tempsample]
    adata_ATACtemp=adata_ATAC[adata_ATAC.obs['sample']==tempsample]
    TG_pseudobulk_temp,RE_pseudobulk_temp=pseudo_bulk(adata_RNAtemp,adata_ATACtemp,singlepseudobulk)                
    TG_pseudobulk=pd.concat([TG_pseudobulk, TG_pseudobulk_temp], axis=1)
    RE_pseudobulk=pd.concat([RE_pseudobulk, RE_pseudobulk_temp], axis=1)
    RE_pseudobulk[RE_pseudobulk > 100] = 100

if not os.path.exists('data/'):
    os.mkdir('data/')
adata_ATAC.write('data/adata_ATAC.h5ad')
adata_RNA.write('data/adata_RNA.h5ad')
TG_pseudobulk=TG_pseudobulk.fillna(0)
RE_pseudobulk=RE_pseudobulk.fillna(0)
pd.DataFrame(adata_ATAC.var['gene_ids']).to_csv('data/Peaks.txt',header=None,index=None)
TG_pseudobulk.to_csv('data/TG_pseudobulk.tsv')
RE_pseudobulk.to_csv('data/RE_pseudobulk.tsv')



method='LINGER'
Datadir='/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/'# This directory should be the same as Datadir defined in the above 'Download the general gene regulatory network' section
GRNdir=Datadir+'data_bulk/'
genome='hg38'
outdir=Datadir + "output/"#output dir
preprocess(TG_pseudobulk,RE_pseudobulk,GRNdir,genome,method,outdir)

activef='ReLU' # active function chose from 'ReLU','sigmoid','tanh'
LINGER_tr.training(GRNdir,method,outdir,activef,'Human')
LL_net.TF_RE_binding(GRNdir,adata_RNA,adata_ATAC,genome,method,outdir)

LL_net.cis_reg(GRNdir,adata_RNA,adata_ATAC,genome,method,outdir)
LL_net.trans_reg(GRNdir,method,outdir,genome)
celltype='all'

LL_net.cell_type_specific_TF_RE_binding(GRNdir,adata_RNA,adata_ATAC,genome,celltype,outdir,method)# different from the previous version
LL_net.cell_type_specific_cis_reg(GRNdir,adata_RNA,adata_ATAC,genome,celltype,outdir)
LL_net.cell_type_specific_trans_reg(GRNdir,adata_RNA,celltype,outdir)
