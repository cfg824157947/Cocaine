import numpy as np
import time
import pandas as pd
from scipy.sparse import csc_matrix, csr_matrix, find
import psutil
import os
from anndata import AnnData
from scipy import sparse
from harmony import harmonize
import scanpy as sc
import scanpy.external as sce


# Make ref from vijay file

file_name = "/data2/duren_lab/cham/cocain/vijay_raw.h5ad"
#file_name = "/project/zduren/durenlab/palmetto/cham/cocain/Matrix/data/vijay_raw.h5ad"
adata = sc.read_h5ad(file_name)
adata.obs['sample'] = adata.obs_names.str.split("-").str[1]
matching = pd.read_csv("/data2/duren_lab/cham/cocain/ncbi_dataset.tsv", sep='\t')
matching['Locus tag'] = matching['Locus tag'].str.replace('_','-',regex=True)
matching = matching[matching['Locus tag'].isin(adata.var_names)]
adata = adata[:,adata.var_names.isin(matching['Locus tag'])]
matching['Locus tag']=pd.Categorical(matching['Locus tag'], categories = adata.var_names)
matching = matching.sort_values(by='Locus tag')

adata.var['Symbol']= matching['Symbol'].values

#adata.var_names = adata.var_names.str.split("-").str[1]

lookup_table = {
    "0": "Optic lobe neurons",
    "1": "5-HT and DA neurons",
    "2": "Glutamatergic neurons",
    "3": "neuropeptides/cholinergic neurons, central brain",
    "4": "Glutamatergic neurons",
    "5": "surface glia",
    "6": "neuropeptides/cholinergic neurons",
    "7": "optic lobe neurons",
    "8": "GABAergic neurons",
    "9": "GABAergic neurons",
    "10": "olfactory projection neurons",
    "11": "optic lobe and antennal lobe neurons",
    "12": "Kenyon cells",
    "13": "Glutamatergic neurons",
    "14": "GABAergic neurons",
    "15": "eye, cuticle, glia",
    "16": "GABAergic neurons",
    "17": "cholinergic neurons of central brain",
    "18": "undetermined, perhaps optic lobe neurons",
    "19": "undetermined",
    "20": "glutamaterigic neurons of central brain",
    "21": "GABAergic neurons of central brain",
    "22": "undetermined",
    "23": "Astrocytes",
    "24": "Astrocytes",
    "25": "Kenyon cells",
    "26": "GABAergic neurons",
    "27": "undetermined",
    "28": "Cells of ellipsoid body and fan shaped body",
    "29": "Photoreceptor cells",
    "30": "GABAergic neurons",
    "31": "undetermined",
    "32": "Mostly undetermined, dorsal fan-shaped body",
    "33": "Surface glia and fat body",
    "34": "undetermined",
    "35": "ellipsoid body, mushroom body",
    "36": "Dopaminergic neurons",
    "37": "undetermined",
    "38": "neuropeptides in perhaps cholinergic neurons",
    "39": "developmental genes",
    "40": "undetermined",
    "41": "undetermined",
    "42": "Tachykinin/neuropeptidergic neurons"
  }
adata.obs["celltypeC"] = adata.obs["celltype"].astype(str)
cluster_names = adata.obs["celltypeC"]
cell_types = [lookup_table[cluster] for cluster in cluster_names]
adata.obs["cell_types"] = cell_types
R = True
if R :
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    PCN = 15
    sc.tl.pca(adata, n_comps=PCN)
    adata.obsm['X_pca'] = pd.read_csv("/data2/duren_lab/cham/cocain/cocaine_Vijay_PCA.csv", index_col=0).values
    Fpca = pd.read_csv("/data2/duren_lab/cham/cocain/feature_PC.csv", index_col=0)
    Fpca = Fpca[Fpca.index.isin(adata.var_names)]
    var_names = Fpca.index
    adata= adata[:, adata.var_names.isin(var_names)]
    adata.var["highly_variable"] = np.ones(adata.shape[1]).astype(bool)
    adata.varm['PCs'] = Fpca.values 
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")
    sc.tl.umap(adata)
    sc.pl.umap(adata, color="cell_types", save=("cell_types_N="+str(15)+"R.pdf"), legend_fontsize=8, legend_loc='right margin', legend_fontweight='normal') #plotting UMAP
else:    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:,adata.var["highly_variable"]] 
    PCN = 15
    sc.tl.pca(adata, n_comps=PCN)
    sc.pp.neighbors(adata, n_neighbors=15)
    sce.pp.harmony_integrate(adata, 'sample', basis ='X_pca_norm' ,max_iter_harmony = 20)
    adata.obsm['X_pca_norm']=((adata.obsm['X_pca'] - adata.obsm['X_pca'].mean(axis=0)) / adata.obsm['X_pca'].std(axis=0))
    list = [15,20,25,30,35,40,45,50]
    list = [15]
    for i in list:
        sc.pp.neighbors(adata, n_neighbors=i, use_rep="X_pca_norm")
        sc.tl.umap(adata)
        sc.pl.umap(adata, color="cell_types", save=("cell_types_N="+str(i)+".pdf"), legend_fontsize=8, legend_loc='right margin', legend_fontweight='normal') #plotting UMAP

adata_ref = adata.copy()
adata_ref.obs['annotation'] = adata_ref.obs['cell_types']
adata_ref.uns["annotation_colors"] = adata_ref.uns["cell_types_colors"]  # fix colors
adata_ref.var_names = adata_ref.var["Symbol"]

ref_file = '/data2/duren_lab/cham/cocain/adata_ref_vijay.h5ad'
adata_ref.write(ref_file)

# Make New adata file

adata = sc.read_10x_h5("/data2/duren_lab/cham/cocain/filtered_feature_bc_matrix.h5",gex_only=False)
adata.obs['sample'] = adata.obs_names.str.split("-").str[1]
adata_RNA = adata[:,adata.var['feature_types']=='Gene Expression']
adata_RNA.var['mt'] = adata_RNA.var_names.str.startswith('mt:')
sc.pp.calculate_qc_metrics(adata_RNA, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(
    adata_RNA,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.1,
    multi_panel=True,
    save = "_check_rna.pdf"
)

keep=['9','23']
adata_RNA = adata_RNA[(adata_RNA.obs.pct_counts_mt < 6) & 
(adata_RNA.obs.n_genes_by_counts< 1600) & 
(adata_RNA.obs.n_genes_by_counts > 200) & 
(adata_RNA.obs.total_counts < 4000) & 
(~adata_RNA.obs['sample'].isin(keep)),:]



adata = adata[adata.obs_names.isin(adata_RNA.obs_names)]
adata_ATAC = adata[:,adata.var['feature_types']=='Peaks']

sc.pp.calculate_qc_metrics(adata_ATAC, percent_top=None, log1p=False, inplace=True)

adata_ATAC= adata_ATAC[(adata_ATAC.obs.n_genes_by_counts< 1800) & 
    (adata_ATAC.obs.total_counts > 400) & 
    (adata_ATAC.obs.total_counts < 4000) ,:]

sc.pl.violin(
    adata_ATAC,
    ["n_genes_by_counts", "total_counts"],
    jitter=0.1,
    multi_panel=True,
    save = "_check_atac.pdf"
)

adata = adata[adata.obs_names.isin(adata_ATAC.obs_names)]
adata_RNA = adata_RNA[adata_RNA.obs_names.isin(adata_ATAC.obs_names)]

# Normalize

sc.pp.normalize_total(adata_RNA, target_sum=1e4)
sc.pp.log1p(adata_RNA, base=2)

sc.pp.highly_variable_genes(adata_RNA, min_mean=0.0125, max_mean=3, min_disp=0.25)
adata_RNA = adata_RNA[:,adata_RNA.var.highly_variable]
sc.pp.scale(adata_RNA, max_value=10)
sc.tl.pca(adata_RNA, n_comps=15)
import scanpy.external as sce
sce.pp.harmony_integrate(adata_RNA,'sample')
adata_RNA.obsm['X_pca_norm']=((adata_RNA.obsm['X_pca_harmony'] - adata_RNA.obsm['X_pca_harmony'].mean(axis=0)) / adata_RNA.obsm['X_pca_harmony'].std(axis=0))
sc.pp.neighbors(adata_RNA, use_rep='X_pca_norm') 
sc.tl.leiden(adata_RNA, resolution = 1, flavor="igraph")
sc.tl.umap(adata_RNA) #embedding
sc.pl.umap(adata_RNA, color="leiden",save="_rna_rna.pdf")

from muon import atac as ac
ac.pp.binarize(adata_ATAC)
sc.pp.filter_genes(adata_ATAC, min_cells=500)
ac.pp.tfidf(adata_ATAC, scale_factor=1e4)
sc.tl.pca(adata_ATAC, n_comps=15)
sce.pp.harmony_integrate(adata_ATAC,'sample')
adata_ATAC.obsm['X_pca_norm']=((adata_ATAC.obsm['X_pca_harmony'] - adata_ATAC.obsm['X_pca_harmony'].mean(axis=0)) / adata_ATAC.obsm['X_pca_harmony'].std(axis=0))
sc.pp.neighbors(adata_ATAC, use_rep='X_pca_norm') 


sc.tl.leiden(adata_ATAC, resolution = 1, flavor="igraph")
sc.tl.umap(adata_ATAC) #embedding
sc.pl.umap(adata_ATAC, color="leiden",save="_atac_atac.pdf")

adata_ATAC.obs['leiden_RNA']=adata_RNA.obs['leiden']
sc.pl.umap(adata_ATAC, color="leiden_RNA",save="_atac_rna.pdf")

adata_RNA.obs['leiden_ATAC']=adata_ATAC.obs['leiden']
sc.pl.umap(adata_RNA, color="leiden_ATAC",save="_rna_atac.pdf")

H=np.vstack((adata_ATAC.obsm['X_pca_norm'].T,adata_RNA.obsm['X_pca_norm'].T))

adata_RNA.obsm['H_scREG']=H.T
sc.pp.neighbors(adata_RNA, use_rep='H_scREG',key_added='scREG') 
sc.tl.leiden(adata_RNA, resolution = 1,key_added='scREG_leiden',neighbors_key='scREG', flavor="igraph")

sc.pl.umap(adata_RNA, color="scREG_leiden",save="_rna_joint.pdf")

adata_ATAC.obsm['H_scREG']=H.T
sc.pp.neighbors(adata_ATAC, use_rep='H_scREG',key_added='scREG') 
sc.tl.leiden(adata_ATAC, resolution = 1,key_added='scREG_leiden',neighbors_key='scREG', flavor="igraph")
sc.pl.umap(adata_ATAC, color="scREG_leiden",save="_atac_joint.pdf")

adata.obsm['H_scREG']=H.T
sc.pp.neighbors(adata, use_rep='H_scREG',key_added='scREG',metric='cosine') 
sc.tl.leiden(adata, resolution = 1,key_added='scREG_leiden',neighbors_key='scREG', flavor="igraph")
sc.tl.umap(adata,neighbors_key='scREG') #embedding
sc.pl.umap(adata, color="scREG_leiden",save="_joint_joint.pdf")

RNA_file = '/data2/duren_lab/cham/cocain/adata_RNA.h5ad'
ATAC_file = '/data2/duren_lab/cham/cocain/adata_ATAC.h5ad'
whole_file = '/data2/duren_lab/cham/cocain/adata.h5ad'
adata_RNA.write(RNA_file)
adata_ATAC.write(ATAC_file)
adata.write(whole_file)








