#set some figure parameters for nice display inside jupyternotebooks.
import scanpy as sc
import scipy.sparse as sp
import pandas as pd
from LingerGRN.preprocess import *
from LingerGRN.pseudo_bulk import *
import os
import torch
import subprocess
import LingerGRN.LINGER_tr as LINGER_tr
import LingerGRN.LL_net as LL_net



#%matplotlib inline
#sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(5, 5), facecolor='white')
#sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

h5_file =  "/data2/duren_lab/cham/cocain/filtered_feature_bc_matrix.h5"
label_file = "/data2/duren_lab/cham/cocain/LINGER/label.csv"

adata = sc.read_10x_h5(h5_file, gex_only=False)
label=pd.read_csv(label_file,header=0, index_col=0)

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


LINGER_data_dir =  "/data2/duren_lab/cham/cocain/LINGER/data/"

if not os.path.exists(LINGER_data_dir):
    os.mkdir(LINGER_data_dir)
adata_ATAC.write((LINGER_data_dir + 'adata_ATAC.h5ad'))
adata_RNA.write((LINGER_data_dir + 'adata_RNA.h5ad'))
TG_pseudobulk=TG_pseudobulk.fillna(0)
RE_pseudobulk=RE_pseudobulk.fillna(0)
pd.DataFrame(adata_ATAC.var['gene_ids']).to_csv((LINGER_data_dir + 'Peaks.txt'),header=None,index=None)
TG_pseudobulk.to_csv((LINGER_data_dir + 'TG_pseudobulk.tsv'))
RE_pseudobulk.to_csv((LINGER_data_dir + 'RE_pseudobulk.tsv'))


Datadir= '/data2/duren_lab/cham/cocain/LINGER/'# This directory should be the same as Datadir defined in the above 'Download the general gene regulatory network' section
GRNdir=Datadir+'provide_data/'
genome='dm6'
outdir=  "/data2/duren_lab/cham/cocain/LINGER/output/" #output dir

activef='ReLU' 
method='scNN'
LINGER_tr.get_TSS(GRNdir,genome,200000) # Here, 200000 represent the largest distance of regulatory element to the TG. Other distance is supported
LINGER_tr.RE_TG_dis(outdir)


activef='ReLU' # active function chose from 'ReLU','sigmoid','tanh'
genomemap=pd.read_csv(GRNdir+'genome_map_homer.txt',sep='\t')
genomemap.index=genomemap['genome_short']
species=genomemap.loc[genome]['species_ensembl']
LINGER_tr.training(GRNdir,method,outdir,activef,species)

LL_net.TF_RE_binding(GRNdir,adata_RNA,adata_ATAC,genome,method,outdir)
LL_net.cis_reg(GRNdir,adata_RNA,adata_ATAC,genome,method,outdir)
LL_net.trans_reg(GRNdir,method,outdir,genome)
celltype='all'

command='paste data/Peaks.bed data/Peaks.txt > data/region.txt'
subprocess.run(command, shell=True)
import pandas as pd
genome_map=pd.read_csv(GRNdir+'genome_map_homer.txt',sep='\t',header=0)
genome_map.index=genome_map['genome_short']
command='findMotifsGenome.pl data/region.txt '+'dm6'+' ./. -size given -find '+GRNdir+'all_motif_rmdup_'+genome_map.loc[genome]['Motif']+'> '+outdir+'MotifTarget.bed'
subprocess.run(command, shell=True)

LL_net.cell_type_specific_TF_RE_binding(GRNdir,adata_RNA,adata_ATAC,genome,celltype,outdir,method)# different from the previous version
LL_net.cell_type_specific_cis_reg(GRNdir,adata_RNA,adata_ATAC,genome,celltype,outdir,method)
LL_net.cell_type_specific_trans_reg(GRNdir,adata_RNA,celltype,outdir)
