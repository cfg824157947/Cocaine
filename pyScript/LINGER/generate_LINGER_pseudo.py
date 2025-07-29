import numpy as np
import pandas as pd
import random
import scanpy as sc
#from LingerGRN.immupute_dis import immupute_dis
def tfidf(ATAC):
    O = 1 * (ATAC > 0)
    tf1 = O / (np.ones((O.shape[0], 1)) * np.log(1 + np.sum(O, axis=0))[np.newaxis,:])
    idf = np.log(1 + O.shape[1] / (1 + np.sum(O > 0, axis=1)))
    O1 = tf1 * (idf[:, np.newaxis] * np.ones((1, O.shape[1])))
    O1[np.isnan(O1)] = 0
    RE = O1.T
    return RE
def find_neighbors(adata_RNA,adata_ATAC):
    import scanpy as sc
    K = 20
    #sc.tl.pca(adata_RNA, svd_solver="arpack")
    sc.pp.normalize_total(adata_RNA, target_sum=1e4)
    sc.pp.log1p(adata_RNA)
    sc.pp.highly_variable_genes(adata_RNA, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_RNA.raw=adata_RNA
    adata_RNA = adata_RNA[:, adata_RNA.var.highly_variable]
    sc.pp.scale(adata_RNA, max_value=10)
    sc.tl.pca(adata_RNA, n_comps=15,svd_solver="arpack")
    pca_RNA=adata_RNA.obsm['X_pca']
    sc.pp.log1p(adata_ATAC)
    sc.pp.highly_variable_genes(adata_ATAC, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_ATAC.raw=adata_ATAC
    adata_ATAC = adata_ATAC[:, adata_ATAC.var.highly_variable]
    sc.pp.scale(adata_ATAC, max_value=10, zero_center=True)
    sc.tl.pca(adata_ATAC, n_comps=15,svd_solver="arpack")
    pca_ATAC=adata_ATAC.obsm['X_pca']
    pca = np.concatenate((pca_RNA,pca_ATAC), axis=1)
    adata_RNA.obsm['pca']=pca
    adata_ATAC.obsm['pca']=pca
    #sc.pp.neighbors(adata_RNA, n_neighbors=K, n_pcs=30,use_rep='pca')
    return adata_RNA,adata_ATAC


def pseudo_bulk(adata_RNA,adata_ATAC,singlepseudobulk, random_seed):
    K = 20
    #sc.tl.pca(adata_RNA, svd_solver="arpack")
    sc.pp.normalize_total(adata_RNA, target_sum=1e4)
    sc.pp.log1p(adata_RNA)
    sc.pp.filter_genes(adata_RNA, min_cells=3)
    sc.pp.highly_variable_genes(adata_RNA, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_RNA.raw=adata_RNA
    adata_RNA = adata_RNA[:, adata_RNA.var.highly_variable]
    sc.pp.scale(adata_RNA, max_value=10)
    sc.tl.pca(adata_RNA, n_comps=15,svd_solver="arpack")
    pca_RNA=adata_RNA.obsm['X_pca']
    sc.pp.log1p(adata_ATAC)
    sc.pp.filter_genes(adata_ATAC, min_cells=3)
    sc.pp.highly_variable_genes(adata_ATAC, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_ATAC.raw=adata_ATAC
    adata_ATAC = adata_ATAC[:, adata_ATAC.var.highly_variable]
    sc.pp.scale(adata_ATAC, max_value=10, zero_center=True)
    sc.tl.pca(adata_ATAC, n_comps=15,svd_solver="arpack")
    pca_ATAC=adata_ATAC.obsm['X_pca']
    pca = np.concatenate((pca_RNA,pca_ATAC), axis=1)
    adata_RNA.obsm['pca']=pca
    adata_ATAC.obsm['pca']=pca
    sc.pp.neighbors(adata_RNA, n_neighbors=K, n_pcs=30,use_rep='pca')
    connectivities=(adata_RNA.obsp['distances']>0)
    import random
    label=pd.DataFrame(adata_RNA.obs['label'])
    label.columns=['label']
    label.index=adata_RNA.obs['barcode'].tolist()
    #label=label['label'].values
    cluster=list(set(label['label'].values))
    allindex=[]
    np.random.seed(random_seed)  # Set seed for reproducibility
    for i in range(len(cluster)):
        temp=label[label['label']==cluster[i]].index
        N = len(temp) # Total number of elements
        if N>=10:
            m = int(np.floor(np.sqrt(N)))+1  # Number of elements to sample
            if singlepseudobulk>0:
                m=1
            sampled_elements = random.sample(range(N), m)
            temp=temp[sampled_elements]
            allindex=allindex+temp.tolist()
    connectivities=pd.DataFrame(connectivities.toarray(),index=adata_RNA.obs['barcode'].tolist())
    connectivities=connectivities.loc[allindex].values
    A=(connectivities @ adata_RNA.raw.X.toarray())
    TG_filter1=A/(K-1)
    RE_filter1=(connectivities @ adata_ATAC.raw.X.toarray())/(K-1)
    TG_filter1=pd.DataFrame(TG_filter1.T,columns=allindex,index=adata_RNA.raw.var['gene_ids'].tolist())
    RE_filter1=pd.DataFrame(RE_filter1.T,columns=allindex,index=adata_ATAC.raw.var['gene_ids'].tolist())
    return TG_filter1,RE_filter1



#set some figure parameters for nice display inside jupyternotebooks.
import scipy
import pandas as pd
import anndata
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse import  csc_matrix
def get_adata(matrix,features,barcodes,label):
    ### generate the anndata
    matrix.data=matrix.data.astype(np.float32)
    adata=anndata.AnnData(X= csc_matrix(matrix.T))
    adata.var['gene_ids']=features[1].values
    adata.obs['barcode']=barcodes[0].values
    if len(barcodes[0].values[0].split("-"))==2:
        adata.obs['sample'] = [int(string.split("-")[1]) for string in barcodes[0].values]
    else:
        adata.obs['sample'] = 1
    rows_to_select=features[features[2]=='Gene Expression'].index
    adata_RNA = adata[:,rows_to_select]
    rows_to_select=features[features[2]=='Peaks'].index
    adata_ATAC = adata[:,rows_to_select]
### if you have the label (cell type annotation)
    idx=adata_RNA.obs['barcode'].isin(label['barcode_use'].values.ravel())
    adata_RNA=adata_RNA[idx]
    adata_ATAC=adata_ATAC[idx]
    label.index=label['barcode_use'].values.flatten()
    adata_RNA.obs['label']=label.loc[adata_RNA.obs['barcode']]['label'].values
    #barcode_indices = np.where(np.isin(adata_RNA.obs['barcode'].values, label['barcode_use'].values))[0]
    #adata_ATAC = adata_ATAC[barcode_indices, :]
    adata_ATAC.obs['label']=label.loc[adata_ATAC.obs['barcode']]['label'].values
    adata_RNA.var["mt"] = adata_RNA.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
    adata_RNA, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
    adata_RNA = adata_RNA[adata_RNA.obs.pct_counts_mt < 5, :].copy()
    adata_RNA.var.index=adata_RNA.var['gene_ids'].values
    adata_RNA.var_names_make_unique()
    adata_RNA.var['gene_ids']=adata_RNA.var.index
    selected_barcode=list(set(adata_RNA.obs['barcode'].values)&set(adata_ATAC.obs['barcode'].values))
    barcode_idx=pd.DataFrame(range(adata_RNA.shape[0]), index=adata_RNA.obs['barcode'].values)
    adata_RNA = adata_RNA[barcode_idx.loc[selected_barcode][0]]
    barcode_idx=pd.DataFrame(range(adata_ATAC.shape[0]), index=adata_ATAC.obs['barcode'].values)
    adata_ATAC = adata_ATAC[barcode_idx.loc[selected_barcode][0]]
    return adata_RNA,adata_ATAC




h5_file =  "/project/zduren/durenlab/cham/cocain/filtered_feature_bc_matrix.h5"
label_file = "/project/zduren/durenlab/cham/cocain/LINGER/label.csv"

adata = sc.read_10x_h5(h5_file, gex_only=False)
label=pd.read_csv(label_file,header=0, index_col=0)

matrix=adata.X.T
adata.var['gene_ids']=adata.var.index
features=pd.DataFrame(adata.var['gene_ids'].values.tolist(),columns=[1])
features[2]=adata.var['feature_types'].values
barcodes=pd.DataFrame(adata.obs_names,columns=[0])
adata_RNA,adata_ATAC=get_adata(matrix,features,barcodes,label)# adata_RNA and adata_ATAC are scRNA and scat Ac

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

#random_seed = 1
sc_metadata_file = '/project/zduren/durenlab/cham/cocain/Data/meta_data3.tsv'
metadata = pd.read_csv(sc_metadata_file, sep='\t', index_col=0)
pseudo_dir = "/project/zduren/durenlab/palmetto/cham/cocain/Final/PseudoBulk/LINGER/"

for random_seed in range(100):
    TG_pseudobulk = pd.DataFrame()  # ← reset
    RE_pseudobulk = pd.DataFrame()  # ← reset
    for tempsample in samplelist:
        adata_RNAtemp = adata_RNA[adata_RNA.obs['sample'] == tempsample].copy()
        adata_ATACtemp = adata_ATAC[adata_ATAC.obs['sample'] == tempsample].copy()

#        adata_RNAtemp=adata_RNA[adata_RNA.obs['sample']==tempsample]
#        adata_ATACtemp=adata_ATAC[adata_ATAC.obs['sample']==tempsample]
        TG_pseudobulk_temp,RE_pseudobulk_temp=pseudo_bulk(adata_RNAtemp,adata_ATACtemp,singlepseudobulk,random_seed)                
        TG_pseudobulk=pd.concat([TG_pseudobulk, TG_pseudobulk_temp], axis=1)
        RE_pseudobulk=pd.concat([RE_pseudobulk, RE_pseudobulk_temp], axis=1)
        RE_pseudobulk[RE_pseudobulk > 100] = 100
    
    
    TG_pseudobulk=TG_pseudobulk.fillna(0)
    RE_pseudobulk=RE_pseudobulk.fillna(0)

    pseudo_TG_path = f"{pseudo_dir}/{random_seed}_TG_pseudo.csv"
    pseudo_RE_path = f"{pseudo_dir}/{random_seed}_RE_pseudo.csv"
    pseudo_meta_path = f"{pseudo_dir}/{random_seed}_meta_pseudo.csv"

    TG_pseudobulk.to_csv(pseudo_TG_path)
    RE_pseudobulk.to_csv(pseudo_RE_path)
    metadata.to_csv(pseudo_meta_path)
    print(f"Pseudo-bulk data for random seed {random_seed} saved.")

