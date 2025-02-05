import pandas as pd
import numpy as np   
from statsmodels.stats.multitest import multipletests
import os




def assignLabel(W,p):
    W=W/(W.sum(axis=0)+1**(-6))
    W2=W.T/(W.T.sum(axis=0)+ 10**(-4));
    max_values = np.max(W2, axis=0)
    max_indices = np.argmax(W2, axis=0)
    quantile = np.percentile(max_values, p*100)
    S_gene=W[:,0];
    S_gene[:]=0
    for i in range(K):
        S_gene[(max_values>quantile)&(max_indices==i)]=i+1
    return S_gene,W2 


def select_group(meta_data,male,pref):
    sex_idx = meta_data['male'] == male
    Line_idx = meta_data['prefer'] == pref
    selected_index = meta_data[sex_idx & Line_idx].index
    return selected_index

def transform_R_analysis_pseudobulk(R_pseudobulk,celltype_list):
    pseudobulk_dict = {}
    for celltype in celltype_list:
        sufix = celltype + "-"
        celltype_cols = R_pseudobulk.filter(like=sufix)  # Selects columns containing 'Astrocyte-'
        celltype_cols.columns = celltype_cols.columns.str.extract(r'(\d+)$')[0]
        pseudobulk_dict[celltype] = celltype_cols

    return pseudobulk_dict

def calculate_logFC_R_analysis_pseudobulk(pseudobulk_dict,metadata):
    celltype_list = pseudobulk_dict.keys()
    FC_df = pd.DataFrame(columns=celltype_list)
    for celltype in pseudobulk_dict.keys():
        print(celltype)
        print(pseudobulk_dict[celltype])
        tmp = metadata[metadata['celltype']==celltype]
        Heroin_idx = tmp['sample_name'][tmp['group']==1]
        Control_idx = tmp['sample_name'][tmp['group']==0]
        Mean1 = pseudobulk_dict[celltype][Heroin_idx.astype(str)].mean(axis=1)
        Mean0 = pseudobulk_dict[celltype][Control_idx.astype(str)].mean(axis=1)
        FC = pd.DataFrame(Mean1/Mean0)
        FC_df[celltype] = FC[0]
    return FC_df


def calculate_logFC(Exp_TG,aud_idx,celltype):
#    aud_idx1=aud_idx.reset_index()
    celltype_barcode = aud_idx[aud_idx['celltype']==celltype].index
    aud_barcode = aud_idx[aud_idx['group']==1].index
    C_barcode = aud_idx[aud_idx['group']==0].index
    celltype_aud_barcode = celltype_barcode[celltype_barcode.isin(aud_barcode)]
    celltype_ctl_barcode = celltype_barcode[~celltype_barcode.isin(aud_barcode)]
#    celltype_aud_barcode = celltype_barcode[celltype_barcode.isin(aud_barcode)]
#    celltype_ctl_barcode = celltype_barcode[~celltype_barcode.isin(aud_barcode)]
    Mean1=Exp_TG[celltype_aud_barcode].mean(axis=1)+10**(-6)
    Mean0=Exp_TG[celltype_ctl_barcode].mean(axis=1)+10**(-6)
    FC=pd.DataFrame(Mean1-Mean0)
    return FC



def calculate_logFC_old(Exp_TG,aud_idx,celltype):
#    aud_idx1=aud_idx.reset_index()
    celltype_barcode = aud_idx[aud_idx['celltype']==celltype].index
    aud_barcode = aud_idx[aud_idx['group']==1].index
    C_barcode = aud_idx[aud_idx['group']==0].index
    celltype_aud_barcode = celltype_barcode[celltype_barcode.isin(aud_barcode)]
    celltype_ctl_barcode = celltype_barcode[~celltype_barcode.isin(aud_barcode)]
#    celltype_aud_barcode = celltype_barcode[celltype_barcode.isin(aud_barcode)]
#    celltype_ctl_barcode = celltype_barcode[~celltype_barcode.isin(aud_barcode)]
    Mean1=Exp_TG[celltype_aud_barcode].mean(axis=1)+10**(-6)
    Mean0=Exp_TG[celltype_ctl_barcode].mean(axis=1)+10**(-6)
    FC=pd.DataFrame(Mean1-Mean0)
    return FC

def make_logFC_df(Exp_TG,aud_idx,celltype_list):
    FC_df = pd.DataFrame()
    for celltype in celltype_list:
        FC = calculate_logFC(Exp_TG,aud_idx,celltype)
        FC_df[celltype] = FC[0]
    return FC_df


def correlation_FC(x,y,method):
    from scipy import stats
    # Loop through each column of y and calculate correlation with x
    correlations = []
    correlationsp=[]
    for i in range(y.shape[1]):
        if method=='pearsonr':
            r, p = stats.pearsonr(x.ravel(), y.values[:,i])  
        if method=='spearmanr':
            r, p = stats.spearmanr(x.ravel(), y.values[:,i]) 
        correlations.append(r)
        correlationsp.append(p)
    correlations=pd.DataFrame(correlations,index=y.columns)
    correlationsp=pd.DataFrame(correlationsp,index=y.columns)
    return correlations,correlationsp

def make_Gene_score_df(DEG_dir,score):
    Gene_score_df = pd.DataFrame()
    Gene_order=None

    for entry in os.listdir(DEG_dir):
        celltype = entry.split('.')[0]
        DEG_file = os.path.join(DEG_dir, entry)
        ranked_genes = pd.read_csv(DEG_file, index_col=2)  # Replace with your ranked list file
        if Gene_order is None:
            Gene_order = ranked_genes.index
            Gene_score_df.index = Gene_order
        ranked_genes = ranked_genes[[score]]  # Replace with your ranked list file
        Gene_score_df[celltype] = ranked_genes.loc[Gene_order]

    return Gene_score_df

def make_Gene_score_df_Drosophira(DEG_dict,score,compaire):
    Gene_score_df = pd.DataFrame()
    Gene_order=None
    for celltype in DEG_dict.keys():
        ranked_genes = DEG_dict[celltype]
        ranked_genes = ranked_genes[ranked_genes['X'].str.contains(compaire)]
        ranked_genes.index = ranked_genes['Gene']
        if Gene_order is None:
            Gene_order = ranked_genes.index
            Gene_score_df.index = Gene_order
        ranked_genes = ranked_genes[[score]]  # Replace with your ranked list file
        Gene_score_df[celltype] = ranked_genes.loc[Gene_order]
    return Gene_score_df


def driver_score_drosophira(reg,adjust_method,corr_method, Gene_score_dict, compaire, score):

    if reg.shape[1]<4:
        reg = reg.drop_duplicates(subset=['RE', 'TF'])
        reg = reg.pivot(index='RE', columns='TF', values='score').fillna(0)
    reg = reg.fillna(0)

    # normalize reg
    cols=reg.sum(axis=0).values
    rows=reg.sum(axis=1).values
    E=np.reshape(rows,(rows.shape[0],1))*np.reshape(cols,(1,cols.shape[0]))/rows.sum()
    #print(E.mean().mean()*10**(-4))
    E=E+E.mean().mean()*10**(-4)
    reg=(reg-E)/E
    reg[reg<0]=0
    reg = reg.loc[~reg.index.duplicated()]
    C_result = pd.DataFrame(index=reg.columns, columns=Gene_score_dict.keys())
    P_result = pd.DataFrame(index=reg.columns, columns=Gene_score_dict.keys())
    Q_result = pd.DataFrame(index=reg.columns, columns=Gene_score_dict.keys())
    # each cell type have different gene list, so we need to loop for
    for celltype in Gene_score_dict.keys():
        print('cell type '+ celltype)
        Gene_score = Gene_score_dict[celltype]
        Gene_score = Gene_score[Gene_score['X'].str.contains(compaire)]
        Gene_score.index = Gene_score['Gene']
        overlap=list(set(Gene_score.index)&set(reg.index))
#        reg=reg.loc[overlap]
        Gene_score=Gene_score.loc[overlap]
        reg=reg.loc[Gene_score.index]


        print(np.isnan(Gene_score.sum()))
        c,cp=correlation_FC(np.log(Gene_score[score]).fillna(0).values,reg,corr_method)
#        C_result=pd.concat([C_result,c],axis=1)
        C_result[celltype] = c
        P_result[celltype] = cp
#        P_result=pd.concat([P_result,cp],axis=1)   
        cp=cp.fillna(1)  
        adjusted_p_values = pd.DataFrame(multipletests(cp[0].values, method=adjust_method)[1],index=c.index)
        Q_result[celltype] = adjusted_p_values
#        pd.concat([Q_result,adjusted_p_values],axis=1)

    return C_result,P_result,Q_result


def diff_Module_cham_simple_drosophila(metadata,celltype_list,S_TG,K,Gene_score_dict,compaire,score):
    pvalue_all= np.zeros((K, len(celltype_list)))
    pvalue_all= pd.DataFrame(pvalue_all,index=['M'+str(i+1) for i in range(K)],columns=celltype_list)
    tvalue_all= np.zeros((K, len(celltype_list)))
    tvalue_all= pd.DataFrame(tvalue_all,index=['M'+str(i+1) for i in range(K)],columns=celltype_list)
#    for k in range(len(celltype)):
    for celltype in celltype_list:
#        temp=Exp_TG.iloc[:,metadata['celltype'].values==celltype]
 #       aud_idxtemp=metadata[(metadata['celltype'].values==celltype)]['group'].values
#        key='C'+celltype+'_anovas'
        key='C'+str(celltype)+'_anovas'
        Gene_score=Gene_score_dict[key]
        Gene_score = Gene_score[Gene_score['X'].str.contains(compaire)]
        Gene_score.index = Gene_score['Gene']
        for i in range(1,K+1):
            print("Module",i)
            Module_genes = S_TG[S_TG['Module']==i]
            overlap=list(set(Gene_score.index)&set(Module_genes.index))
            adj_p_values = Gene_score.loc[overlap, score]
            pvalue_all.loc['M'+str(i),celltype]= (adj_p_values < 0.1).sum() / len(adj_p_values)
            #            adj_p_values = DEG_score_dict[celltype].loc[S_TG['Module']==i]
#            adj_p_values = DEG_score_dict[celltype].loc[S_TG['Module']==i, celltype]
        # Print the results
    return pvalue_all


def driver_score_cham(reg, Gene_score, adjust_method,corr_method):

    if reg.shape[1]<4:
        reg = reg.drop_duplicates(subset=['RE', 'TF'])
        reg = reg.pivot(index='RE', columns='TF', values='score').fillna(0)
    reg = reg.fillna(0)

    # normalize reg
    cols=reg.sum(axis=0).values
    rows=reg.sum(axis=1).values
    E=np.reshape(rows,(rows.shape[0],1))*np.reshape(cols,(1,cols.shape[0]))/rows.sum()
    #print(E.mean().mean()*10**(-4))
    E=E+E.mean().mean()*10**(-4)
    reg=(reg-E)/E
    reg[reg<0]=0
    reg = reg.loc[~reg.index.duplicated()]
    overlap=list(set(Gene_score.index)&set(reg.index))
    reg=reg.loc[overlap]
    Gene_score=Gene_score.loc[overlap]
    reg=reg.loc[overlap]
    reg=reg.loc[Gene_score.index]

    C_result = pd.DataFrame(index=reg.columns, columns=Gene_score.columns)
    P_result = pd.DataFrame(index=reg.columns, columns=Gene_score.columns)
    Q_result = pd.DataFrame(index=reg.columns, columns=Gene_score.columns)
    for celltype in Gene_score.columns:
 #       print('cell type '+ celltype)
#        Gene_score=Gene_score.loc[reg.index]
#        print(np.isnan(Gene_score.sum()))
#        c,cp=correlation_FC(np.log(Gene_score[celltype]).fillna(0).values,reg,corr_method)
        c,cp=correlation_FC((Gene_score[celltype]).fillna(0).values,reg,corr_method)
#        C_result=pd.concat([C_result,c],axis=1)
        C_result[celltype] = c
        P_result[celltype] = cp
#        P_result=pd.concat([P_result,cp],axis=1)   
        cp=cp.fillna(1)  
        adjusted_p_values = pd.DataFrame(multipletests(cp[0].values, method=adjust_method)[1],index=c.index)
        Q_result[celltype] = adjusted_p_values
#        pd.concat([Q_result,adjusted_p_values],axis=1)

    return C_result,P_result,Q_result


def driver_result(C_result,Q_result,K):
    import pandas as pd
    # Create a sample DataFrame of size 100x7   
    # Rank all values in the DataFrame
    top_5_rows=[]
    for column in C_result.columns:
        sorted_column = C_result[column].sort_values(ascending=False)
        top_5_rows =top_5_rows+sorted_column.index[:K].tolist()
    top_5_rows=list(set(top_5_rows))
    for column in C_result.columns:
        sorted_column = C_result[column].sort_values(ascending=True)
        top_5_rows =top_5_rows+sorted_column.index[:K].tolist()
    top_5_rows=list(set(top_5_rows))
    #df_ranked = np.abs(C_result).max(axis=1).rank(ascending=False)
    # Get the indices of the top 20 values
    #top_20_indices = df_ranked.sort_values()[:K].index
    # Extract the row and column indices from the multi-index
    #top_20_rows = top_20_indices.copy()
    #cutoff=np.abs(C_result.loc[top_20_rows]).max(axis=1).min()
    #top_20_rows=C_result.index[(np.abs(C_result)>cutoff).sum(axis=1)>0]
    return C_result.loc[top_5_rows],Q_result.loc[top_5_rows]


def assignLabel(W,p,K):
    W=W/(W.sum(axis=0)+1**(-6))
    W2=W.T/(W.T.sum(axis=0)+ 10**(-4));
    max_values = np.max(W2, axis=0)
    max_indices = np.argmax(W2, axis=0)
    quantile = np.percentile(max_values, p*100)
    S_gene=W[:,0];
    S_gene[:]=0
    for i in range(K):
        S_gene[(max_values>quantile)&(max_indices==i)]=i+1
    return S_gene,W2 

def diff_Module_cham(metadata,celltype_list,S_TG,K,DEG_df):
    pvalue_all= np.zeros((K, len(celltype_list)))
    pvalue_all= pd.DataFrame(pvalue_all,index=['M'+str(i+1) for i in range(K)],columns=celltype_list)
    tvalue_all= np.zeros((K, len(celltype_list)))
    tvalue_all= pd.DataFrame(tvalue_all,index=['M'+str(i+1) for i in range(K)],columns=celltype_list)
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    from scipy.stats import ttest_1samp
#    for k in range(len(celltype)):
    for celltype in celltype_list:
#        temp=Exp_TG.iloc[:,metadata['celltype'].values==celltype]
 #       aud_idxtemp=metadata[(metadata['celltype'].values==celltype)]['group'].values
        for i in range(1,K+1):
            print("Module",i)
            logfc_values = DEG_df.loc[S_TG['Module']==i, celltype]
            stat, p_value = ttest_1samp(logfc_values, 0)
            tvalue_all.loc['M'+str(i),celltype]=stat
            pvalue_all.loc['M'+str(i),celltype]=p_value
        # Print the results
    return pvalue_all,tvalue_all


def diff_Module_cham_simple(metadata,celltype_list,S_TG,K,DEG_df):
    pvalue_all= np.zeros((K, len(celltype_list)))
    pvalue_all= pd.DataFrame(pvalue_all,index=['M'+str(i+1) for i in range(K)],columns=celltype_list)
    tvalue_all= np.zeros((K, len(celltype_list)))
    tvalue_all= pd.DataFrame(tvalue_all,index=['M'+str(i+1) for i in range(K)],columns=celltype_list)
#    for k in range(len(celltype)):
    for celltype in celltype_list:
#        temp=Exp_TG.iloc[:,metadata['celltype'].values==celltype]
 #       aud_idxtemp=metadata[(metadata['celltype'].values==celltype)]['group'].values
        for i in range(1,K+1):
            print("Module",i)
            adj_p_values = DEG_df.loc[S_TG['Module']==i, celltype]
            pvalue_all.loc['M'+str(i),celltype]= (adj_p_values < 0.1).sum() / len(adj_p_values)
        # Print the results
    return pvalue_all



def Module_trans_cham(trans_reg,metadata,Gene_score_df,K_list,p_list,save_dir,GWASfile=None,GWAS_score=None,simple=False):
    import numpy as np
    from scipy import stats


    overlap=list(set(Gene_score_df.index)&set(trans_reg.index))
    Gene_score_df=Gene_score_df.loc[overlap]
    trans_reg=trans_reg.loc[Gene_score_df.index]
    TFset=trans_reg.columns
    TGset=trans_reg.index

    R1 = stats.zscore(trans_reg,1);R1[np.isnan(R1)] = 0.0
    R2 = stats.zscore(trans_reg,0);R2[np.isnan(R2)] = 0.0
    from sklearn.preprocessing import quantile_transform
    Z=R1+R2
    Z[Z<0]=0  
    if GWAS_score is not None:
        alpha = 1
        TF_GWAS_score = pd.DataFrame((~Z.columns.isin(GWAS_score.index)).astype(int), index = Z.columns)
        tmp_Gidx = TF_GWAS_score[TF_GWAS_score[0]==0].index
        TF_GWAS_score[TF_GWAS_score[0]==0] = GWAS_score.loc[tmp_Gidx]
        TF_GWAS_score = (1 - (alpha * (1 + np.log10(TF_GWAS_score))))

        TG_GWAS_score = pd.DataFrame((~Z.index.isin(GWAS_score.index)).astype(int), index = Z.index)
        tmp_Gidx = TG_GWAS_score[TG_GWAS_score[0]==0].index
        TG_GWAS_score[TG_GWAS_score[0]==0] = GWAS_score.loc[tmp_Gidx]
        TG_GWAS_score = (1 - (alpha * (1 + np.log10(TG_GWAS_score))))
        #TG_GWAS_score = (1 -np.log10(TG_GWAS_score))/2

        TF_score_matrix = np.diag(TF_GWAS_score.values.flatten())
        TG_score_matrix = np.diag(TG_GWAS_score.values.flatten())

        old_Z = Z
        Z = TG_score_matrix @ Z# @ TF_score_matrix


    celltype_list = metadata['celltype'].unique()
    print('identify modules......')
    from sklearn.decomposition import NMF

    for K in K_list:
        nmf = NMF(n_components=K, init='random', random_state=0,max_iter=1000)
        W = nmf.fit_transform(Z)
        H = nmf.components_
        for p in p_list:
            [S_TG,W2]=assignLabel(W,p,K);
            [S_TF,H2]=assignLabel(H.T,p,K);
            S_TG=pd.DataFrame(S_TG,index=TGset,columns=['Module'])
            S_TF=pd.DataFrame(S_TF,index=TFset,columns=['Module'])
            
            print('differential modules......')
            if simple:
                significant_num = diff_Module_cham_simple(metadata,celltype_list,S_TG,K,Gene_score_df)
                nlog_P = significant_num
                p_value_file =save_dir + "/Module_genes/K"+ str(K) +"p" + str(p) +"pvalue.csv"
            else:
                pvalue_all,tvalue_all=diff_Module_cham(metadata,celltype_list,S_TG,K,Gene_score_df)
                nlog_P = -np.log10(pvalue_all)
                p_value_file =save_dir + "/Module_genes/K"+ str(K) +"p" + str(p) +"nlog10_p.csv"

            if not os.path.exists(save_dir+ "/Module_genes/"):
                os.makedirs(save_dir+ "/Module_genes/")

            nlog_P.to_csv(p_value_file)
            S_TG.to_csv(save_dir + "/Module_genes/K"+ str(K) +"p" + str(p) +"TG.csv")
            S_TF.to_csv(save_dir +  "/Module_genes/K"+ str(K) +"p" + str(p) +"TF.csv")

            figure_dir = save_dir + '/figure/LINGER_module_compare/'
            if not os.path.exists(figure_dir):
                os.makedirs(figure_dir)

            plt.figure(figsize=(10, 6))
            sns.heatmap(nlog_P.T, cmap="viridis")
            plt.title("K"+ str(K) +"p" + str(p))
            plt.savefig((figure_dir + ("K"+ str(K) +"p" + str(p)) + "Module_genes.png"),dpi=300)
            plt.clf()




def Module_trans_cham_drosophira(trans_reg,metadata,Gene_score_dict,compaire,score,K_list,p_list,save_dir,GWASfile=None,GWAS_score=None,simple=True):
    import numpy as np
    from scipy import stats


    TFset=trans_reg.columns
    TGset=trans_reg.index

    R1 = stats.zscore(trans_reg,1);R1[np.isnan(R1)] = 0.0
    R2 = stats.zscore(trans_reg,0);R2[np.isnan(R2)] = 0.0
    from sklearn.preprocessing import quantile_transform
    Z=R1+R2
    Z[Z<0]=0  
    if GWAS_score is not None:
        alpha = 1
        TF_GWAS_score = pd.DataFrame((~Z.columns.isin(GWAS_score.index)).astype(int), index = Z.columns)
        tmp_Gidx = TF_GWAS_score[TF_GWAS_score[0]==0].index
        TF_GWAS_score[TF_GWAS_score[0]==0] = GWAS_score.loc[tmp_Gidx]
        TF_GWAS_score = (1 - (alpha * (1 + np.log10(TF_GWAS_score))))

        TG_GWAS_score = pd.DataFrame((~Z.index.isin(GWAS_score.index)).astype(int), index = Z.index)
        tmp_Gidx = TG_GWAS_score[TG_GWAS_score[0]==0].index
        TG_GWAS_score[TG_GWAS_score[0]==0] = GWAS_score.loc[tmp_Gidx]
        TG_GWAS_score = (1 - (alpha * (1 + np.log10(TG_GWAS_score))))
        #TG_GWAS_score = (1 -np.log10(TG_GWAS_score))/2

        TF_score_matrix = np.diag(TF_GWAS_score.values.flatten())
        TG_score_matrix = np.diag(TG_GWAS_score.values.flatten())

        old_Z = Z
        Z = TG_score_matrix @ Z# @ TF_score_matrix


    celltype_list = metadata['celltype'].unique()
    celltype_list = metadata['cluster'].unique()
    print('identify modules......')
    from sklearn.decomposition import NMF

    for K in K_list:
        nmf = NMF(n_components=K, init='random', random_state=0,max_iter=1000)
        W = nmf.fit_transform(Z)
        H = nmf.components_
        for p in p_list:
            [S_TG,W2]=assignLabel(W,p,K);
            [S_TF,H2]=assignLabel(H.T,p,K);
            S_TG=pd.DataFrame(S_TG,index=TGset,columns=['Module'])
            S_TF=pd.DataFrame(S_TF,index=TFset,columns=['Module'])
            
            print('differential modules......')
            if simple:
                significant_num = diff_Module_cham_simple_drosophila(metadata,celltype_list,S_TG,K,Gene_score_dict,compaire,score)

                nlog_P = significant_num
                p_value_file =save_dir + "/Module_genes/K"+ str(K) +"p" + str(p) +"pvalue.csv"
                #p_value_file =save_dir + "/Module_genes/"+ compaire+ "_" +score + "_K"+ str(K) +"p" + str(p) +"pvalue.csv"
            else:
                pvalue_all,tvalue_all=diff_Module_cham(metadata,celltype_list,S_TG,K,Gene_score_df)
                nlog_P = -np.log10(pvalue_all)
                p_value_file =save_dir + "/Module_genes/K"+ str(K) +"p" + str(p) +"nlog10_p.csv"

            if not os.path.exists(save_dir+ "/Module_genes/"):
                os.makedirs(save_dir+ "/Module_genes/")

            nlog_P.to_csv(p_value_file)
            S_TG.to_csv(save_dir + "/Module_genes/K"+ str(K) +"p" + str(p) +"TG.csv")
            S_TF.to_csv(save_dir +  "/Module_genes/K"+ str(K) +"p" + str(p) +"TF.csv")

            figure_dir = save_dir + '/figure/LINGER_module_compare/'
            if not os.path.exists(figure_dir):
                os.makedirs(figure_dir)

            plt.figure(figsize=(10, 6))
            sns.heatmap(nlog_P.T, cmap="viridis")
            plt.title("K"+ str(K) +"p" + str(p))
            plt.savefig((figure_dir + ("K"+ str(K) +"p" + str(p)) + "Module_genes.png"),dpi=300)
            plt.clf()


def driver_score_cham(reg, Gene_score, adjust_method='bonferroni',corr_method='pearsonr'):
    if reg.shape[1]<4:
        reg = reg.drop_duplicates(subset=['RE', 'TF'])
        reg = reg.pivot(index='RE', columns='TF', values='score').fillna(0)
    reg = reg.fillna(0)
    # normalize reg
    cols=reg.sum(axis=0).values
    rows=reg.sum(axis=1).values
    E=np.reshape(rows,(rows.shape[0],1))*np.reshape(cols,(1,cols.shape[0]))/rows.sum()
    #print(E.mean().mean()*10**(-4))
    E=E+E.mean().mean()*10**(-4)
    reg=(reg-E)/E
    reg[reg<0]=0
    reg = reg.loc[~reg.index.duplicated()]
    overlap=list(set(Gene_score.index)&set(reg.index))
    Gene_score=Gene_score.loc[overlap]
    reg=reg.loc[Gene_score.index]

    C_result = pd.DataFrame(index=reg.columns, columns=Gene_score.columns)
    P_result = pd.DataFrame(index=reg.columns, columns=Gene_score.columns)
    Q_result = pd.DataFrame(index=reg.columns, columns=Gene_score.columns)
    for celltype in Gene_score.columns:
#        print('cell type '+ celltype)
#        Gene_score=Gene_score.loc[reg.index]
#        print(np.isnan(Gene_score.sum()))
        c,cp=correlation_FC(Gene_score[celltype].fillna(0).values,reg,corr_method)
        C_result[celltype] = c
        P_result[celltype] = cp
        cp=cp.fillna(1)  
        adjusted_p_values = pd.DataFrame(multipletests(cp[0].values, method=adjust_method)[1],index=c.index)
        Q_result[celltype] = adjusted_p_values

    return C_result,P_result,Q_result
