#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import scipy.sparse as sp
import scipy.io
import anndata as ad
import gzip
import shutil

def convert_tsv_to_adata_and_sparse(tsv_path, h5ad_path, outdir, feature_type="peak"):
    df = pd.read_csv(tsv_path, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'barcode', 'count'])
    df['region'] = df['chrom'].astype(str) + ":" + df['start'].astype(str) + "-" + df['end'].astype(str)
    df = df[['region', 'barcode', 'count']]

    row_index = pd.Categorical(df['region'])
    col_index = pd.Categorical(df['barcode'])

    X = sp.coo_matrix(
        (df['count'], (row_index.codes, col_index.codes)),
        shape=(len(row_index.categories), len(col_index.categories))
    ).tocsr()

    adata = ad.AnnData(X.T)
    adata.var_names = row_index.categories
    adata.obs_names = col_index.categories
    adata.write(h5ad_path)

    os.makedirs(outdir, exist_ok=True)
    scipy.io.mmwrite(f"{outdir}/matrix.mtx", adata.X.T)
    pd.Series(adata.obs_names).to_csv(f"{outdir}/barcodes.tsv", index=False, header=False)

    features = pd.DataFrame({
        "feature_id": adata.var_names,
        "feature_name": adata.var_names,
        "feature_type": feature_type
    })
    features.to_csv(f"{outdir}/features.tsv", sep="\t", index=False, header=False)

    for fn in ["matrix.mtx", "barcodes.tsv", "features.tsv"]:
        with open(f"{outdir}/{fn}", 'rb') as f_in:
            with gzip.open(f"{outdir}/{fn}.gz", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(f"{outdir}/{fn}")

    print(f"✅ Done: {h5ad_path} and sparse folder: {outdir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert BED-style peak count TSV to h5ad and sparse matrix")
    parser.add_argument('--tsv', required=True, help="Input .tsv file with chrom, start, end, barcode, count")
    parser.add_argument('--h5ad', required=True, help="Output path for .h5ad")
    parser.add_argument('--outdir', required=True, help="Output directory for sparse matrix")
    parser.add_argument('--feature_type', default="peak", help="Feature type for features.tsv")

    args = parser.parse_args()
    convert_tsv_to_adata_and_sparse(args.tsv, args.h5ad, args.outdir, args.feature_type)
