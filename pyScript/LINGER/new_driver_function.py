import os
import numpy as np
import pandas as pd

from statsmodels.stats.multitest import multipletests
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import NMF


def assignLabel(W, p):
    """
    Assign module labels for either genes (TG) or regulators (TF) using
    NMF component matrix W. A label is assigned if the maximum
    normalized value across components is above a quantile threshold.

    Parameters
    ----------
    W : np.ndarray
        A 2D matrix (e.g., from NMF component).
    p : float
        Quantile to apply as a threshold (0 < p < 1).
    K : int
        Number of components (modules).
    
    Returns
    -------
    S : np.ndarray
        1D array of assigned module labels. 0 indicates no assignment,
        1..K indicate module labels.
    W2 : np.ndarray
        The row-normalized version of W for each entity.
    """
    # Column normalize W
    W = W / (W.sum(axis=0) + 1e-6)
    K = W.shape[1]
    # Row normalize
    W2 = (W.T / (W.T.sum(axis=0) + 1e-4)).T

    # Find the max component per row (i.e., which module it belongs to)
    max_values = np.max(W2, axis=1)
    max_indices = np.argmax(W2, axis=1)

    # Determine threshold
    quantile = np.percentile(max_values, p * 100)

    # Create label vector
    S = np.zeros(W2.shape[0], dtype=int)
    for i in range(K):
        idx = (max_values > quantile) & (max_indices == i)
        S[idx] = i + 1

    return S, W2


def select_group(meta_data, male, pref):
    """
    Subset the meta_data DataFrame by 'male' and 'prefer' columns.

    Parameters
    ----------
    meta_data : pd.DataFrame
        Metadata containing columns ['male', 'prefer'].
    male : int/bool
        Value to match in 'male' column.
    pref : str
        Value to match in 'prefer' column.

    Returns
    -------
    selected_index : pd.Index
        The index of rows matching the given male/pref criteria.
    """
    sex_idx = meta_data['male'] == male
    Line_idx = meta_data['prefer'] == pref
    selected_index = meta_data[sex_idx & Line_idx].index
    return selected_index

def select_group_by_key(meta_data, key, value):
    """
    Subset the meta_data DataFrame by a key-value pair.

    Parameters
    ----------
    meta_data : pd.DataFrame
        Metadata containing
    """
    idx = meta_data[key] == value
    selected_index = meta_data[idx].index
    return selected_index

def select_group(meta_data,male=None,pref=None,Line=None):
    """
    Subset the meta_data DataFrame by 'male' and 'prefer' columns.

    Parameters
    ----------
    meta_data : pd.DataFrame
        Metadata containing columns ['male', 'prefer'].
    male : int/bool
        Value to match in 'male' column.
    pref : str
        Value to match in 'prefer' column.
    Line: str
        Value to match in 'Line' column.

    Returns
    -------
    selected_index : pd.Index
        The index of rows matching the given male/pref criteria.
    """
    if male is not None:
        sex_idx = meta_data['male'] == male
    if pref is not None:
        Line_idx = meta_data['prefer'] == pref
    if Line is not None:
        Line_idx = meta_data['Line'] == Line
    selected_index = meta_data[sex_idx & Line_idx].index
    return selected_index

def transform_R_analysis_pseudobulk(R_pseudobulk, celltype_list):
    """
    Subset a pseudobulk DataFrame for each celltype, renaming columns
    by sample name only.

    Parameters
    ----------
    R_pseudobulk : pd.DataFrame
        Pseudobulk expression matrix with column naming convention like 'Astrocyte-10'.
    celltype_list : list of str
        List of celltypes to subset.

    Returns
    -------
    pseudobulk_dict : dict
        Dictionary with key=celltype, value=subset DataFrame of columns
        for that celltype.
    """
    pseudobulk_dict = {}
    for celltype in celltype_list:
        suffix = celltype + "-"
        celltype_cols = R_pseudobulk.filter(like=suffix)  # e.g. 'Astrocyte-'
        # Extract only the numeric part of the column name
        celltype_cols.columns = celltype_cols.columns.str.extract(r'(\d+)$')[0]
        pseudobulk_dict[celltype] = celltype_cols
    return pseudobulk_dict


def calculate_logFC_R_analysis_pseudobulk(pseudobulk_dict, metadata):
    """
    Calculate fold-change per gene in each celltype using
    the pseudobulk approach.

    Parameters
    ----------
    pseudobulk_dict : dict
        Dictionary where keys=celltype, values=expression DataFrame.
    metadata : pd.DataFrame
        Must contain ['celltype', 'sample_name', 'group'] columns, where
        'group' has 1 and 0 for two conditions.

    Returns
    -------
    FC_df : pd.DataFrame
        Fold-change for each gene across all celltypes.
    """
    celltype_list = pseudobulk_dict.keys()
    FC_df = pd.DataFrame(columns=celltype_list)
    for celltype in celltype_list:
        tmp = metadata[metadata['celltype'] == celltype]
        heroin_idx = tmp['sample_name'][tmp['group'] == 1].astype(str)
        control_idx = tmp['sample_name'][tmp['group'] == 0].astype(str)

        mean1 = pseudobulk_dict[celltype][heroin_idx].mean(axis=1)
        mean0 = pseudobulk_dict[celltype][control_idx].mean(axis=1)

        FC_df[celltype] = mean1 / mean0  # ratio
    return FC_df


def calculate_logFC(exp_df, meta_df, celltype):
    """
    Compute log fold change (difference) between group=1 and group=0
    for a given celltype.

    Parameters
    ----------
    exp_df : pd.DataFrame
        Expression matrix (rows=genes, cols=cells).
    meta_df : pd.DataFrame
        Cell-level metadata with columns ['celltype', 'group'] at least.
    celltype : str
        The celltype of interest.

    Returns
    -------
    FC : pd.DataFrame
        A single-column DataFrame with the difference in mean expression
        for each gene in the given celltype.
    """
    celltype_barcode = meta_df[meta_df['celltype'] == celltype].index
    aud_barcode = meta_df[meta_df['group'] == 1].index
    ctl_barcode = meta_df[meta_df['group'] == 0].index

    celltype_aud_barcode = celltype_barcode.intersection(aud_barcode)
    celltype_ctl_barcode = celltype_barcode.intersection(ctl_barcode)

    mean1 = exp_df[celltype_aud_barcode].mean(axis=1) + 1e-6
    mean0 = exp_df[celltype_ctl_barcode].mean(axis=1) + 1e-6
    FC = pd.DataFrame(mean1 - mean0, columns=['logFC'])
    return FC


def make_logFC_df(exp_df, meta_df, celltype_list):
    """
    Create a DataFrame of log fold changes (differences) for all genes
    in each celltype.

    Parameters
    ----------
    exp_df : pd.DataFrame
        Expression matrix (rows=genes, cols=cells).
    meta_df : pd.DataFrame
        Cell-level metadata with columns ['celltype', 'group'] at least.
    celltype_list : list
        List of celltypes to compute logFC.

    Returns
    -------
    FC_df : pd.DataFrame
        Columns are celltypes, rows are genes. Each entry is the difference
        in mean expression for group=1 vs group=0.
    """
    FC_df = pd.DataFrame()
    for celltype in celltype_list:
        FC = calculate_logFC(exp_df, meta_df, celltype)
        FC_df[celltype] = FC['logFC']
    return FC_df


def correlation_FC(x, y, method):
    """
    Compute correlation of a single vector x with each column in a DataFrame y.

    Parameters
    ----------
    x : np.array or list
        1D array to correlate against each column of y.
    y : pd.DataFrame
        Each column in y will be correlated with x.
    method : str
        'pearsonr' or 'spearmanr'.

    Returns
    -------
    correlations : pd.DataFrame
        r-values for each column in y.
    correlationsp : pd.DataFrame
        p-values for each column in y.
    """
    from scipy import stats

    correlations = []
    correlationsp = []
    for i in range(y.shape[1]):
        col_data = y.iloc[:, i].values
        if method == 'pearsonr':
            r, p = stats.pearsonr(x, col_data)
        elif method == 'spearmanr':
            r, p = stats.spearmanr(x, col_data)
        else:
            raise ValueError("Unsupported correlation method: choose 'pearsonr' or 'spearmanr'")
        correlations.append(r)
        correlationsp.append(p)

    correlations = pd.DataFrame(correlations, index=y.columns)
    correlationsp = pd.DataFrame(correlationsp, index=y.columns)
    return correlations, correlationsp


def make_Gene_score_df(DEG_dir, score):
    """
    Create a DataFrame of gene scores from multiple CSV files in a directory.
    Each CSV is assumed to correspond to a particular celltype.

    Parameters
    ----------
    DEG_dir : str
        Path to the directory containing CSV files.
    score : str
        Column name to extract from each CSV.

    Returns
    -------
    Gene_score_df : pd.DataFrame
        Rows=genes, columns=celltypes, entries=score.
    """
    Gene_score_df = pd.DataFrame()
    Gene_order = None

    for entry in os.listdir(DEG_dir):
        if not entry.endswith('.csv'):
            continue
        celltype = entry.split('.')[0]
        DEG_file = os.path.join(DEG_dir, entry)
        ranked_genes = pd.read_csv(DEG_file, index_col=2)
        if Gene_order is None:
            Gene_order = ranked_genes.index
            Gene_score_df.index = Gene_order
        # Grab the requested score column
        sub_scores = ranked_genes[[score]]  # single-column DataFrame
        Gene_score_df[celltype] = sub_scores.loc[Gene_order, score]

    return Gene_score_df


def make_Gene_score_df_Drosophira(DEG_dict, score, compaire):
    """
    Similar to make_Gene_score_df but for Drosophila, using a dictionary
    of DataFrames (DEG_dict) and subsetting by 'compaire' in column 'X'.

    Parameters
    ----------
    DEG_dict : dict
        Keys=celltypes, values=DataFrames of differential expression.
    score : str
        Column name to use as the score.
    compaire : str
        String pattern to subset rows from each celltype's DataFrame.

    Returns
    -------
    Gene_score_df : pd.DataFrame
        Rows=genes, columns=celltypes, entries=score.
    """
    Gene_score_df = pd.DataFrame()
    Gene_order = None

    for celltype, df in DEG_dict.items():
        ranked_genes = df[df['X'].str.contains(compaire)]
        ranked_genes.index = ranked_genes['Gene']

        if Gene_order is None:
            Gene_order = ranked_genes.index
            Gene_score_df.index = Gene_order

        sub_scores = ranked_genes[[score]]
        Gene_score_df[celltype] = sub_scores.loc[Gene_order, score]

    return Gene_score_df


def driver_score_drosophira(reg, adjust_method, corr_method, Gene_score_dict, compaire, score):
    """
    Compute correlation between a regulatory matrix reg and gene scores
    for multiple celltypes in a dictionary (Drosophila version).

    Parameters
    ----------
    reg : pd.DataFrame
        Regulatory scores. Expected columns=TFs, index=REs (genes).
    adjust_method : str
        Method passed to statsmodels.stats.multitest.multipletests.
    corr_method : str
        'pearsonr' or 'spearmanr'.
    Gene_score_dict : dict
        Keys=celltypes, values=DataFrame with columns=['X', 'Gene', <score>].
    compaire : str
        Pattern used to filter rows in each celltype DataFrame by 'X'.
    score : str
        Column name in the gene score DataFrame to correlate.

    Returns
    -------
    C_result, P_result, Q_result : pd.DataFrame
        Correlation coefficients, raw p-values, and adjusted p-values
        (across TFs in columns, celltypes in columns).
    """
    if reg.shape[1] < 4:
        # pivot if the matrix is in 'RE, TF, score' format
        reg = reg.drop_duplicates(subset=['RE', 'TF'])
        reg = reg.pivot(index='RE', columns='TF', values='score').fillna(0)
    reg = reg.fillna(0)

    # Normalize reg
    cols = reg.sum(axis=0).values
    rows = reg.sum(axis=1).values
    E = np.outer(rows, cols) / rows.sum()
    E = E + E.mean() * 1e-4
    reg = (reg - E) / E
    reg[reg < 0] = 0
    reg = reg.loc[~reg.index.duplicated()]

    C_result = pd.DataFrame(index=reg.columns, columns=Gene_score_dict.keys())
    P_result = pd.DataFrame(index=reg.columns, columns=Gene_score_dict.keys())
    Q_result = pd.DataFrame(index=reg.columns, columns=Gene_score_dict.keys())

    for celltype, gene_score_df in Gene_score_dict.items():
        print(f"Processing celltype {celltype}")
        # Filter by 'compaire' in 'X'
        gene_score = gene_score_df[gene_score_df['X'].str.contains(compaire)]
        gene_score.index = gene_score['Gene']

        overlap = set(gene_score.index).intersection(reg.index)
        if not overlap:
            # No overlapping genes, skip
            continue

        # subset
        gene_score_sub = gene_score.loc[overlap, score].fillna(0)
        reg_sub = reg.loc[overlap]

        c, cp = correlation_FC(np.log(gene_score_sub + 1e-6).values, reg_sub, corr_method)

        C_result[celltype] = c.squeeze()
        P_result[celltype] = cp.squeeze()

        cp = cp.fillna(1)
        adjusted_p = multipletests(cp[0].values, method=adjust_method)[1]
        Q_result[celltype] = pd.DataFrame(adjusted_p, index=c.index)[0]

    return C_result, P_result, Q_result


def diff_Module_cham_simple_drosophila(metadata, celltype_list, S_TG, K, Gene_score_dict, compaire, score):
    """
    Simple proportion-based differential module test for Drosophila.
    For each module, compute proportion of genes with a certain threshold
    (e.g., adj p-value < 0.1) in the given gene score dictionary.

    Parameters
    ----------
    metadata : pd.DataFrame
        Cell-level metadata (not heavily used here).
    celltype_list : list
        All celltypes to consider.
    S_TG : pd.DataFrame
        Module assignments for genes: columns=['Module'].
    K : int
        Number of modules.
    Gene_score_dict : dict
        Keys are strings like 'C<celltype>_anovas', values are DataFrames.
    compaire : str
        String pattern to subset rows by 'X'.
    score : str
        Column name in each DataFrame to interpret as p-value or score.

    Returns
    -------
    pvalue_all : pd.DataFrame
        Row=Module, col=celltype. Each entry is the proportion of genes
        within that module that pass the threshold (score < 0.1).
    """
    pvalue_all = pd.DataFrame(np.zeros((K, len(celltype_list))),
                              index=[f'M{i+1}' for i in range(K)],
                              columns=celltype_list)

    for celltype in celltype_list:
        key = 'C' + str(celltype) + '_anovas'
        gene_score = Gene_score_dict[key]
        gene_score = gene_score[gene_score['X'].str.contains(compaire)]
        gene_score.index = gene_score['Gene']

        for i in range(1, K+1):
            module_genes = S_TG[S_TG['Module'] == i].index
            overlap = module_genes.intersection(gene_score.index)
            if len(overlap) == 0:
                continue
            values = gene_score.loc[overlap, score]
            # proportion with value < 0.1
            pvalue_all.loc[f'M{i}', celltype] = np.mean(values < 0.1)

    return pvalue_all


def driver_score_cham(reg, Gene_score, adjust_method='bonferroni', corr_method='pearsonr'):
    """
    Correlate a regulatory matrix (reg) with a gene score matrix (Gene_score).
    Both should share the same row index (genes). Each column in Gene_score
    is treated as a separate condition/celltype.

    Parameters
    ----------
    reg : pd.DataFrame
        Regulatory scores, shape=(genes, TFs). If shape[1]<4,
        it is assumed to be melted form with columns=['RE','TF','score'].
    Gene_score : pd.DataFrame
        shape=(genes, celltypes). Each celltype is a column of gene-level scores.
    adjust_method : str
        Method for multiple testing correction. E.g. 'bonferroni', 'fdr_bh'.
    corr_method : str
        'pearsonr' or 'spearmanr'.

    Returns
    -------
    C_result : pd.DataFrame
        Correlation coefficients. Index=TFs, columns=celltypes.
    P_result : pd.DataFrame
        Raw p-values.
    Q_result : pd.DataFrame
        Adjusted p-values with the specified method.
    """
    # If needed, pivot from melted to wide
    if reg.shape[1] < 4:
        reg = reg.drop_duplicates(subset=['RE', 'TF'])
        reg = reg.pivot(index='RE', columns='TF', values='score').fillna(0)

    reg = reg.fillna(0)

    # Normalize
    cols = reg.sum(axis=0).values
    rows = reg.sum(axis=1).values
    E = np.outer(rows, cols) / rows.sum()
    E = E + E.mean() * 1e-4
    reg = (reg - E) / E
    reg[reg < 0] = 0
    reg = reg.loc[~reg.index.duplicated()]

    # Overlap
    overlap = list(set(Gene_score.index).intersection(reg.index))
    Gene_score = Gene_score.loc[overlap]
    reg = reg.loc[overlap]

    # Allocate results
    C_result = pd.DataFrame(index=reg.columns, columns=Gene_score.columns)
    P_result = pd.DataFrame(index=reg.columns, columns=Gene_score.columns)
    Q_result = pd.DataFrame(index=reg.columns, columns=Gene_score.columns)

    # Correlate each celltype in Gene_score vs. reg
    for celltype in Gene_score.columns:
        x = Gene_score[celltype].fillna(0).values
        c, cp = correlation_FC(x, reg, corr_method)

        C_result[celltype] = c.squeeze()
        P_result[celltype] = cp.squeeze()

        cp = cp.fillna(1)
        adjusted_p = multipletests(cp[0].values, method=adjust_method)[1]
        Q_result[celltype] = pd.DataFrame(adjusted_p, index=c.index)[0]

    return C_result, P_result, Q_result


def diff_Module_cham(metadata, celltype_list, S_TG, K, DEG_df):
    """
    Perform a simple t-test (per module) for whether gene-level logFC
    is different from zero.

    Parameters
    ----------
    metadata : pd.DataFrame
        Not heavily used, but provides celltypes.
    celltype_list : list
        All celltypes to consider.
    S_TG : pd.DataFrame
        Module assignments for genes: columns=['Module'].
    K : int
        Number of modules.
    DEG_df : pd.DataFrame
        Per-gene logFC per celltype (rows=genes, columns=celltypes).

    Returns
    -------
    pvalue_all : pd.DataFrame
        p-values from one-sample t-tests across modules.
    tvalue_all : pd.DataFrame
        t-statistics from the same test.
    """
    from scipy.stats import ttest_1samp

    pvalue_all = pd.DataFrame(np.zeros((K, len(celltype_list))),
                              index=[f'M{i+1}' for i in range(K)],
                              columns=celltype_list)

    tvalue_all = pd.DataFrame(np.zeros((K, len(celltype_list))),
                              index=[f'M{i+1}' for i in range(K)],
                              columns=celltype_list)

    for celltype in celltype_list:
        for i in range(1, K + 1):
            module_genes = S_TG.index[S_TG['Module'] == i]
            logfc_values = DEG_df.loc[module_genes, celltype].dropna()
            if len(logfc_values) == 0:
                continue

            stat, p_value = ttest_1samp(logfc_values, 0)
            tvalue_all.loc[f'M{i}', celltype] = stat
            pvalue_all.loc[f'M{i}', celltype] = p_value

    return pvalue_all, tvalue_all


def diff_Module_cham_simple(metadata, celltype_list, S_TG, K, DEG_df, threshold=0.1):
    """
    Compute proportion of genes in each module that pass a significance
    threshold in their corresponding columns (e.g., p-value < 0.1).

    Parameters
    ----------
    metadata : pd.DataFrame
        Cell-level metadata (not heavily used).
    celltype_list : list
        All celltypes to consider.
    S_TG : pd.DataFrame
        Module assignments for genes: columns=['Module'].
    K : int
        Number of modules.
    DEG_df : pd.DataFrame
        Per-gene p-value or adjusted p-value (rows=genes, columns=celltypes).
    threshold : float
        Default 0.1. The cutoff for counting “significant” or “passed”.

    Returns
    -------
    pvalue_all : pd.DataFrame
        Proportions. Rows=modules, columns=celltypes.
    """
    pvalue_all = pd.DataFrame(np.zeros((K, len(celltype_list))),
                              index=[f'M{i+1}' for i in range(K)],
                              columns=celltype_list)

    for celltype in celltype_list:
        for i in range(1, K + 1):
            module_genes = S_TG.index[S_TG['Module'] == i]
            values = DEG_df.loc[module_genes, celltype].dropna()

            if len(values) == 0:
                continue

            # The proportion of genes that have a score < threshold
            pvalue_all.loc[f'M{i}', celltype] = np.mean(values < threshold)

    return pvalue_all


def Module_trans_cham(trans_reg, metadata, Gene_score_df, K_list, p_list,
                      save_dir, GWAS_score=None, simple=False):
    """
    Identify modules (via NMF on a z-scored + clipped version of trans_reg),
    assign modules to TGs and TFs, and then run differential module analysis.

    Parameters
    ----------
    trans_reg : pd.DataFrame
        Regulatory matrix (rows=genes, columns=TFs).
    metadata : pd.DataFrame
        Cell-level or sample-level metadata with 'celltype' column.
    Gene_score_df : pd.DataFrame
        Per-gene scores, typically logFC or p-value, shape=(genes, celltypes).
    K_list : list
        List of integers for the number of modules (NMF components).
    p_list : list
        List of floats for module assignment quantile thresholds.
    save_dir : str
        Output directory to store results.
    GWAS_score : pd.Series, optional
        If provided, modifies the matrix by weighting rows & columns.
    simple : bool
        If True, uses proportion-based test (diff_Module_cham_simple).
        If False, uses t-tests (diff_Module_cham).
    """
    from sklearn.preprocessing import quantile_transform

    # Overlap
    overlap = set(Gene_score_df.index).intersection(trans_reg.index)
    Gene_score_df = Gene_score_df.loc[overlap]
    trans_reg = trans_reg.loc[overlap]

    # Z-score rows/cols
    R1 = stats.zscore(trans_reg, axis=1)
    R1[np.isnan(R1)] = 0.0

    R2 = stats.zscore(trans_reg, axis=0)
    R2[np.isnan(R2)] = 0.0

    Z = R1 + R2
    Z[Z < 0] = 0

    # Optional weighting by GWAS
    if GWAS_score is not None:
        alpha = 1.0
        # TF side
        tf_gwas_mask = (~Z.columns.isin(GWAS_score.index)).astype(int)
        TF_GWAS_score = pd.DataFrame(tf_gwas_mask.values, index=Z.columns)
        # Fill known with actual GWAS
        known_tf_idx = TF_GWAS_score[TF_GWAS_score[0] == 0].index
        TF_GWAS_score.loc[known_tf_idx, 0] = GWAS_score.loc[known_tf_idx]
        TF_GWAS_score = (1 - (alpha * (1 + np.log10(TF_GWAS_score))))

        # TG side
        tg_gwas_mask = (~Z.index.isin(GWAS_score.index)).astype(int)
        TG_GWAS_score = pd.DataFrame(tg_gwas_mask.values, index=Z.index)
        known_tg_idx = TG_GWAS_score[TG_GWAS_score[0] == 0].index
        TG_GWAS_score.loc[known_tg_idx, 0] = GWAS_score.loc[known_tg_idx]
        TG_GWAS_score = (1 - (alpha * (1 + np.log10(TG_GWAS_score))))

        # Build diagonal weighting
        TF_score_matrix = np.diag(TF_GWAS_score[0].values)
        TG_score_matrix = np.diag(TG_GWAS_score[0].values)

        # Weighted version
        Z = TG_score_matrix @ Z  # Optionally could multiply TF side too

    celltype_list = metadata['celltype'].unique()

    if not os.path.exists(os.path.join(save_dir, "Module_genes")):
        os.makedirs(os.path.join(save_dir, "Module_genes"))

    for K in K_list:
        nmf_model = NMF(n_components=K, init='random', random_state=0, max_iter=1000)
        W = nmf_model.fit_transform(Z)
        H = nmf_model.components_

        for p in p_list:
            # Gene modules
            S_TG, _ = assignLabel(W, p, K)
            # TF modules
            S_TF, _ = assignLabel(H.T, p, K)

            S_TG_df = pd.DataFrame({'Module': S_TG}, index=Z.index)
            S_TF_df = pd.DataFrame({'Module': S_TF}, index=Z.columns)

            # Evaluate differential modules
            if simple:
                significant_num = diff_Module_cham_simple(
                    metadata, celltype_list, S_TG_df, K, Gene_score_df
                )
                nlog_P = significant_num
                p_value_file = os.path.join(
                    save_dir, "Module_genes", f"K{K}p{p}proportion.csv"
                )
            else:
                pvalue_all, tvalue_all = diff_Module_cham(
                    metadata, celltype_list, S_TG_df, K, Gene_score_df
                )
                nlog_P = -np.log10(pvalue_all)
                p_value_file = os.path.join(
                    save_dir, "Module_genes", f"K{K}p{p}nlog10_p.csv"
                )

            nlog_P.to_csv(p_value_file)
            S_TG_df.to_csv(os.path.join(save_dir, "Module_genes", f"K{K}p{p}TG.csv"))
            S_TF_df.to_csv(os.path.join(save_dir, "Module_genes", f"K{K}p{p}TF.csv"))

            # Plot
            figure_dir = os.path.join(save_dir, 'figure', 'LINGER_module_compare')
            if not os.path.exists(figure_dir):
                os.makedirs(figure_dir)

            plt.figure(figsize=(10, 6))
            sns.heatmap(nlog_P.T, cmap="viridis")
            plt.title(f"K={K}, p={p}")
            out_png = os.path.join(figure_dir, f"K{K}p{p}_Module_genes.png")
            plt.savefig(out_png, dpi=300)
            plt.clf()


def Module_trans_cham_drosophira(trans_reg, metadata, Gene_score_dict, compaire, score,
                                 K_list, p_list, save_dir, GWAS_score=None, simple=True):
    """
    Similar pipeline for Drosophila data, where Gene_score_dict is a dictionary
    of DataFrames. NMF is performed, modules assigned, then differential modules
    computed by proportion-based approach (default).

    Parameters
    ----------
    trans_reg : pd.DataFrame
        Regulatory matrix (rows=genes, columns=TFs).
    metadata : pd.DataFrame
        Contains celltype/cluster info in 'celltype' or 'cluster'.
    Gene_score_dict : dict
        Keys=some tag, values=DEG-like DataFrames with 'X', 'Gene', <score>.
    compaire : str
        Pattern used to filter each DataFrame in Gene_score_dict by 'X'.
    score : str
        Column name used as the metric or p-value.
    K_list : list
        List of integers (# of modules for NMF).
    p_list : list
        List of floats (quantile thresholds).
    save_dir : str
        Output directory.
    GWAS_score : pd.Series, optional
        GWAS weighting factor for each gene (and possibly each TF).
    simple : bool
        If True, uses diff_Module_cham_simple_drosophila to compute proportions.
    """
    from sklearn.preprocessing import quantile_transform

    TFset = trans_reg.columns
    TGset = trans_reg.index

    # Z-score
    R1 = stats.zscore(trans_reg, axis=1)
    R1[np.isnan(R1)] = 0.0

    R2 = stats.zscore(trans_reg, axis=0)
    R2[np.isnan(R2)] = 0.0

    Z = R1 + R2
    Z[Z < 0] = 0

    # Optionally apply GWAS weighting
    if GWAS_score is not None:
        alpha = 1.0

        tf_gwas_mask = (~Z.columns.isin(GWAS_score.index)).astype(int)
        TF_GWAS_score = pd.DataFrame(tf_gwas_mask.values, index=Z.columns)
        known_tf_idx = TF_GWAS_score[TF_GWAS_score[0] == 0].index
        TF_GWAS_score.loc[known_tf_idx, 0] = GWAS_score.loc[known_tf_idx]
        TF_GWAS_score = (1 - (alpha * (1 + np.log10(TF_GWAS_score))))

        tg_gwas_mask = (~Z.index.isin(GWAS_score.index)).astype(int)
        TG_GWAS_score = pd.DataFrame(tg_gwas_mask.values, index=Z.index)
        known_tg_idx = TG_GWAS_score[TG_GWAS_score[0] == 0].index
        TG_GWAS_score.loc[known_tg_idx, 0] = GWAS_score.loc[known_tg_idx]
        TG_GWAS_score = (1 - (alpha * (1 + np.log10(TG_GWAS_score))))

        TF_score_matrix = np.diag(TF_GWAS_score[0].values)
        TG_score_matrix = np.diag(TG_GWAS_score[0].values)
        Z = TG_score_matrix @ Z

    # celltype_list
    celltype_list = metadata['celltype'].unique() if 'celltype' in metadata.columns \
        else metadata['cluster'].unique()

    if not os.path.exists(os.path.join(save_dir, "Module_genes")):
        os.makedirs(os.path.join(save_dir, "Module_genes"))

    for K in K_list:
        nmf_model = NMF(n_components=K, init='random', random_state=0, max_iter=1000)
        W = nmf_model.fit_transform(Z)
        H = nmf_model.components_

        for p in p_list:
            S_TG, _ = assignLabel(W, p, K)
            S_TF, _ = assignLabel(H.T, p, K)

            S_TG_df = pd.DataFrame({'Module': S_TG}, index=TGset)
            S_TF_df = pd.DataFrame({'Module': S_TF}, index=TFset)

            # Proportion-based differential modules (Drosophila)
            if simple:
                significant_num = diff_Module_cham_simple_drosophila(
                    metadata, celltype_list, S_TG_df, K, Gene_score_dict, compaire, score
                )
                nlog_P = significant_num
                p_value_file = os.path.join(
                    save_dir, "Module_genes", f"K{K}p{p}_pvalue.csv"
                )
            else:
                # Or a different approach if you want (not implemented in this example)
                raise NotImplementedError("Set simple=True for Drosophila approach.")

            nlog_P.to_csv(p_value_file)
            S_TG_df.to_csv(os.path.join(save_dir, "Module_genes", f"K{K}p{p}TG.csv"))
            S_TF_df.to_csv(os.path.join(save_dir, "Module_genes", f"K{K}p{p}TF.csv"))

            figure_dir = os.path.join(save_dir, 'figure', 'LINGER_module_compare')
            if not os.path.exists(figure_dir):
                os.makedirs(figure_dir)

            plt.figure(figsize=(10, 6))
            sns.heatmap(nlog_P.T, cmap="viridis")
            plt.title(f"Drosophila K={K}, p={p}")
            out_png = os.path.join(figure_dir, f"K{K}p{p}_Module_genes.png")
            plt.savefig(out_png, dpi=300)
            plt.clf()

def net_col_permute(network):
    """
    Permute columns of a network DataFrame.
    """
    #net work is a pandas dataframe
    Permu_net = network.copy()
    for i in range(Permu_net.shape[1]):
        Permu_net.iloc[:,i] = np.random.permutation(Permu_net.iloc[:,i])
    return Permu_net

def driver_score_cham_new(reg, Gene_score, adjust_method='bonferroni', corr_method='pearsonr'):
    """
    Correlate a regulatory matrix (reg) with a gene score matrix (Gene_score).
    Both should share the same row index (genes). Each column in Gene_score
    is treated as a separate condition/celltype.

    Parameters
    ----------
    reg : pd.DataFrame
        Regulatory scores, shape=(genes, TFs). If shape[1]<4,
        it is assumed to be melted form with columns=['RE','TF','score'].
    Gene_score : pd.DataFrame
        shape=(genes, celltypes). Each celltype is a column of gene-level scores.
    adjust_method : str
        Method for multiple testing correction. E.g. 'bonferroni', 'fdr_bh'.
    corr_method : str
        'pearsonr' or 'spearmanr'.

    Returns
    -------
    C_result : pd.DataFrame
        Correlation coefficients. Index=TFs, columns=celltypes.
    P_result : pd.DataFrame
        Raw p-values.
    Q_result : pd.DataFrame
        Adjusted p-values with the specified method.
    """
    # If needed, pivot from melted to wide
    if reg.shape[1] < 4:
        reg = reg.drop_duplicates(subset=['RE', 'TF'])
        reg = reg.pivot(index='RE', columns='TF', values='score').fillna(0)

    reg = reg.fillna(0)

    # Normalize
    cols = reg.sum(axis=0).values
    rows = reg.sum(axis=1).values
    E = np.outer(rows, cols) / rows.sum()
    E = E + E.mean() * 1e-4
    reg = (reg - E) / E
    reg[reg < 0] = 0
    reg = reg.loc[~reg.index.duplicated()]

    # Overlap
    overlap = list(set(Gene_score.index).intersection(reg.index))
    Gene_score = Gene_score.loc[overlap]
    reg = reg.loc[overlap]

    # Allocate results
    C_result_dict = {'original': pd.DataFrame(0,index=reg.columns, columns=Gene_score.columns),
                     'mean': pd.DataFrame(0,index=reg.columns, columns=Gene_score.columns),
                     'std': pd.DataFrame(0,index=reg.columns, columns=Gene_score.columns),
                     'zscore': pd.DataFrame(0,index=reg.columns, columns=Gene_score.columns)}
    P_result_dict = {'original': pd.DataFrame(0,index=reg.columns, columns=Gene_score.columns),
                     'mean': pd.DataFrame(0,index=reg.columns, columns=Gene_score.columns),
                     'std': pd.DataFrame(0,index=reg.columns, columns=Gene_score.columns),
                     'zscore': pd.DataFrame(0,index=reg.columns, columns=Gene_score.columns)}
    Q_result_dict = {'original': pd.DataFrame(0,index=reg.columns, columns=Gene_score.columns),
                        'mean': pd.DataFrame(0,index=reg.columns, columns=Gene_score.columns),
                        'std': pd.DataFrame(0,index=reg.columns, columns=Gene_score.columns),
                        'zscore': pd.DataFrame(0,index=reg.columns, columns=Gene_score.columns)}



    # caluculate each score with permutation network
    permute = 30
    for celltype in Gene_score.columns:
        x = Gene_score[celltype].fillna(0).values
        for i in range(permute):
            permu_reg = net_col_permute(reg)
            c, cp = correlation_FC(x, permu_reg, corr_method)
            C_result_dict['mean'][celltype] += c.squeeze()
            C_result_dict['std'][celltype] += c.squeeze()**2
            P_result_dict['mean'][celltype] += cp.squeeze()
            P_result_dict['std'][celltype] += cp.squeeze()**2
            cp = cp.fillna(1)
            adjusted_p = multipletests(cp[0].values, method=adjust_method)[1]
            Q_result_dict['mean'][celltype] += pd.DataFrame(adjusted_p, index=c.index)[0]
            Q_result_dict['std'][celltype] += pd.DataFrame(adjusted_p, index=c.index)[0]**2

        C_result_dict['mean'][celltype] /= permute
        C_result_dict['std'][celltype] = np.sqrt(C_result_dict['std'][celltype]/permute - C_result_dict['mean'][celltype]**2)
        P_result_dict['mean'][celltype] /= permute
        P_result_dict['std'][celltype] = np.sqrt(P_result_dict['std'][celltype]/permute - P_result_dict['mean'][celltype]**2)
        Q_result_dict['mean'][celltype] /= permute
        Q_result_dict['std'][celltype] = np.sqrt(Q_result_dict['std'][celltype]/permute - Q_result_dict['mean'][celltype]**2)

        #original network
        c, cp = correlation_FC(x, reg, corr_method)
        C_result_dict['original'][celltype] = c.squeeze()
        P_result_dict['original'][celltype] = cp.squeeze()
        cp = cp.fillna(1)
        adjusted_p = multipletests(cp[0].values, method=adjust_method)[1]
        Q_result_dict['original'][celltype] = pd.DataFrame(adjusted_p, index=c.index)[0]

        #calculate zscore
        C_result_dict['zscore'][celltype] = (C_result_dict['original'][celltype] - C_result_dict['mean'][celltype])/C_result_dict['std'][celltype]
        P_result_dict['zscore'][celltype] = (P_result_dict['original'][celltype] - P_result_dict['mean'][celltype])/P_result_dict['std'][celltype]
        Q_result_dict['zscore'][celltype] = (Q_result_dict['original'][celltype] - Q_result_dict['mean'][celltype])/Q_result_dict['std'][celltype]

    return C_result_dict, P_result_dict, Q_result_dict
