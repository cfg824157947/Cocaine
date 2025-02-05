import LingerGRN.LL_net as LL_net
import LingerGRN.LINGER_tr as LINGER_tr
from LingerGRN.preprocess import *
from LingerGRN.pseudo_bulk import *
import os
import scipy.sparse as sp
import pandas as pd
import scanpy as sc

from LingerGRN.preprocess import *









Datadir='/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/'# This directory should be the same as Datadir defined in the above 'Download the general gene regulatory network' section



adata_ATAC = sc.read_h5ad((Datadir + 'data/adata_ATAC.h5ad'))
adata_RNA = sc.read_h5ad((Datadir + 'data/adata_RNA.h5ad'))
TG_pseudobulk=pd.read_csv((Datadir + 'data/TG_pseudobulk.tsv'),index_col=0)
RE_pseudobulk=pd.read_csv((Datadir + 'data/RE_pseudobulk.tsv'),index_col=0)


#TG=pd.read_csv((Datadir + 'data/TG_pseudobulk.tsv'),index_col=0)

method='LINGER'
GRNdir=Datadir+'data_bulk/'
genome='hg38'
outdir=Datadir + "output/"#output dir
activef='ReLU' # active function chose from 'ReLU','sigmoid','tanh'
celltype='all'



preprocess(TG_pseudobulk,RE_pseudobulk,GRNdir,genome,method,outdir)
print("preprocess Done")
LINGER_tr.training(GRNdir,method,outdir,activef,'Human')
print("tr.training Done")
LL_net.TF_RE_binding(GRNdir,adata_RNA,adata_ATAC,genome,method,outdir)
print("TF_RE_binding Done")

LL_net.cis_reg(GRNdir,adata_RNA,adata_ATAC,genome,method,outdir)
print("cis_reg Done")
LL_net.trans_reg(GRNdir,method,outdir,genome)
print("trans_reg Done")

print("celltype")
LL_net.cell_type_specific_TF_RE_binding(GRNdir,adata_RNA,adata_ATAC,genome,celltype,outdir,method)# different from the previous version
print("TF_RE_binding Done")
LL_net.cell_type_specific_cis_reg(GRNdir,adata_RNA,adata_ATAC,genome,celltype,outdir,method)
print("cis_reg Done")
LL_net.cell_type_specific_trans_reg(GRNdir,adata_RNA,celltype,outdir)
print("trans_reg Done")
