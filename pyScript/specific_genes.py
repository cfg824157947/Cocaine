import numpy as np
import pandas as pd
from anndata import AnnData
import scanpy as sc

xlsx_name = "/data2/duren_lab/cham/cocain/DEG_EXP.xlsx"
adata = sc.read_10x_h5("/data2/duren_lab/cham/cocain/filtered_feature_bc_matrix.h5")
adata_RNA_wholeGene = sc.read_h5ad("/data2/duren_lab/cham/cocain/adata_RNA_wholeGene.h5ad")
metadata = pd.read_csv("./Data/split_meta_data.tsv", sep = '\t', index_col=0)
idx = adata.obs_names.isin(metadata.index)
adata_counts = adata[idx]

adata_log1p.obs['celltype'] = adata_log1p.obs['celltype'].astype('category')
sc.tl.rank_genes_groups(adata_log1p, "celltype", method="t-test")
pval_adj_df = pd.DataFrame(adata_log1p.uns["rank_genes_groups"]["pvals_adj"])
names_df = pd.DataFrame(adata_log1p.uns["rank_genes_groups"]["names"])
logfoldchanges_df = pd.DataFrame(adata_log1p.uns["rank_genes_groups"]["logfoldchanges"])
df_list = []

with pd.ExcelWriter(xlsx_name) as writer:
    for i in range(47):
        parti_adata = adata_counts[adata_counts.obs['celltype']==i]
        gene_means = pd.DataFrame(np.mean(parti_adata.X, axis=0).T, index = parti_adata.var_names)
        df = pd.DataFrame( {'Gene': names_df[str(i)] ,'AvgExp': gene_means.loc[names_df[str(i)].values][0].values, 'pval_adj':pval_adj_df[str(i)].values, 'logFC': logfoldchanges_df[str(i)].values},index = names_df[str(i)])
        df = df.sort_values('pval_adj')
        df.to_excel(writer, sheet_name=(('cluster' + str(i))), index=False)
