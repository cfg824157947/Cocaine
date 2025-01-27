import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
Datadir='/data2/duren_lab/cham/cocain/LINGER/'
outdir=Datadir + "output/"#output dir
metadata = pd.read_csv((Datadir + "meta_data2.csv"), index_col=0)
#metadata['group'] = (metadata['group'] == 'HEROIN').astype(int)
TG_pseudobulk = pd.read_csv((Datadir + "data/TG_pseudobulk.tsv"), index_col=0)
#Gene_score = pd.read_csv("/project/zduren/durenlab/palmetto/cham/Heroin/ori_data/GENE_score_uniq.tsv" ,sep = '\t', index_col=0)
K=10
p=0.5

figure_dir = Datadir+'figure/'
plt.figure(figsize=(30, 12))


def diff_Module(Exp_TG,metadata,celltype,S_TG,K):
    pvalue_all= np.zeros((K, len(celltype)))
    tvalue_all= np.zeros((K, len(celltype)))
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    for k in range(len(celltype)):
        temp=Exp_TG.iloc[:,metadata['celltype'].values==celltype[k]]
        aud_idxtemp=metadata[(metadata['celltype'].values==celltype[k])]['group'].values
        Exp_mean=stats.zscore(temp.T).T.groupby(S_TG['Module'].values).mean()
        Exp_mean=Exp_mean.loc[range(1,K+1)]
        print(celltype[k])
        print(Exp_mean.shape)
        #if Exp_mean.shape[1]<2:
        #    continue
        X=Exp_mean.values[:,(aud_idxtemp==1)]
        Y=Exp_mean.values[:,(aud_idxtemp==0)]
        p_values = np.zeros((K, ))   
        t_values = np.zeros((K, ))
        from scipy.stats import ttest_ind
        for i in range(K):
            print("Module",i)
            t_values[i], p_values[i] = ttest_ind(X[i], Y[i])
        pvalue_all[:,k]=p_values
        tvalue_all[:,k]=t_values
    
    adjusted_p_values = multipletests(p_values, method='fdr_bh')[1]
    pvalue_all=pd.DataFrame(pvalue_all,index=['M'+str(i+1) for i in range(K)],columns=celltype)
    tvalue_all=pd.DataFrame(tvalue_all,index=['M'+str(i+1) for i in range(K)],columns=celltype)
    return pvalue_all,tvalue_all


def Module_trans(outdir,metadata,TG_pseudobulk,K,GWASfile=None):
    import numpy as np
    from scipy import stats

    
    print('loading GRN......')
    
    trans_reg=pd.read_csv(outdir+'cell_population_trans_regulatory.txt',sep='\t',index_col=0)
    #TG_pseudobulk=TG_pseudobulk_all
    TFset=trans_reg.columns
    TGset=trans_reg.index
    R1 = stats.zscore(trans_reg,1);R1[np.isnan(R1)] = 0.0
    R2 = stats.zscore(trans_reg,0);R2[np.isnan(R2)] = 0.0
    #R1[R1<0]=0;R2[R2<0]=0
    #X=TG_pseudobulk.loc[TGset]
    from sklearn.preprocessing import quantile_transform
    #Exp=quantile_transform(np.log2(X + 1), n_quantiles=10, random_state=0, copy=True)
    #Exp=quantile_transform(np.log2(TG_pseudobulk + 1), n_quantiles=10, random_state=0, copy=True)
    Exp=pd.DataFrame(TG_pseudobulk,index=TG_pseudobulk.index,columns=TG_pseudobulk.columns)
    Z=R1+R2
    Z[Z<0]=0  
    print('identify modules......')
    from sklearn.decomposition import NMF
    nmf = NMF(n_components=K, init='random', random_state=0)
    W = nmf.fit_transform(Z)
    H = nmf.components_
    [S_TG,W2]=assignLabel(W,p);
    [S_TF,H2]=assignLabel(H.T,p);
    #Exp_mean=stats.zscore(Exp_TG.T).T.groupby(S_TG).mean()
    S_TG=pd.DataFrame(S_TG,index=TGset,columns=['Module'])
    S_TF=pd.DataFrame(S_TF,index=TFset,columns=['Module'])
    Exp_TF=Exp.loc[TFset]
    Exp_TG=Exp.loc[TGset]
    for m in range(1,K+1):
        S_TG[S_TG['Module'] == m].to_csv(Datadir + "/Module_Genes/K"+ str(K) +"p" + str(p) +"M"+ str(m)  + ".csv")

    print('differential modules......')
    Line_list = metadata['Line'].unique()
    celltype = metadata['celltype'].unique()
    prefer_list = metadata['prefer'].unique()
#    for Line in Line_list:

    for Line in prefer_list:
        for male in [0,1]:
            sex_idx = (metadata['male']==male).values

            if male==1:
                sex='Male'
            else:
                sex='Female'
    #        Line_idx = (metadata['Line']==Line).values
            Line_idx = (metadata['prefer']==Line).values
            metadata_tmp = metadata[Line_idx&sex_idx]
            using_celltype = []
            for cell in celltype:
                groupN = metadata_tmp[metadata_tmp['celltype']==cell]['group'].unique().shape[0]
                if groupN == 2:
                    print(cell)
                    using_celltype = using_celltype + [cell]
            print(Line + "###################################")
            pvalue_all, tvalue_all=diff_Module(Exp_TG.loc[:,Line_idx&sex_idx],metadata[Line_idx&sex_idx],using_celltype,S_TG,K)
            nlog_P = -np.log10(pvalue_all)
            sns.heatmap(nlog_P.T,cmap="viridis")
            plt.title(Line)
            plt.savefig(figure_dir + "/Module_pval_heatmap_"+ sex + "_" + Line + "_K" + str(K) +"_p"+ str(p) +".png",dpi=300)
            plt.clf()



    if GWASfile is not None:
        print('GWAS enrich......')
        GWASgene=pd.DataFrame(index=TG_pseudobulk.index)
        for i in range(len(GWASfile)):
            temp=pd.read_csv(GWASfile[i],sep='\t',header=None)
            idx=np.zeros((GWASgene.shape[0],1))
            idx[GWASgene.index.isin(temp[0].values),:]=1
            GWASgene['GWAS_'+(str(i+1))]=idx
        p_fisher,odds_ratio_fisher=GWAS_Module_enrich(S_TG,TGset,GWASgene,K)
        Module_result=Module_obj()
        Module_result.S_TG=S_TG
        Module_result.pvalue_all=pvalue_all
        Module_result.tvalue_all=tvalue_all
        Module_result.p_fisher=p_fisher
        Module_result.odds_ratio_fisher=odds_ratio_fisher
        return Module_result
    else:
        Module_result=Module_obj()
        Module_result.S_TG=S_TG
        Module_result.pvalue_all=pvalue_all
        Module_result.tvalue_all=tvalue_all
        return Module_result


#def detect_module(trans_reg,metadata,Gene_score_df,K):

def Module_trans_cham(trans_reg,metadata,Gene_score_df,K_list,p_list,save_dir,GWASfile=None,GWAS_score=None,simple=False,surfix=''):
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
                p_value_file =save_dir + "/Module_genes/" + surfix + "K"+ str(K) +"p" + str(p) +"pvalue.csv"
            else:
                pvalue_all,tvalue_all=diff_Module_cham(metadata,celltype_list,S_TG,K,Gene_score_df)
                nlog_P = -np.log10(pvalue_all)
                p_value_file =save_dir + "/Module_genes/" + surfix + "K"+ str(K) +"p" + str(p) +"nlog10_p.csv"

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
            plt.title(surfix + "K"+ str(K) +"p" + str(p))
            plt.savefig((figure_dir + (surfix + "K"+ str(K) +"p" + str(p)) + "Module_genes.png"),dpi=300)
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
        print('cell type '+ celltype)
#        Gene_score=Gene_score.loc[reg.index]
        print(np.isnan(Gene_score.sum()))
        c,cp=correlation_FC(Gene_score[celltype].fillna(0).values,reg,corr_method)
        C_result[celltype] = c
        P_result[celltype] = cp
        cp=cp.fillna(1)  
        adjusted_p_values = pd.DataFrame(multipletests(cp[0].values, method=adjust_method)[1],index=c.index)
        Q_result[celltype] = adjusted_p_values

    return C_result,P_result,Q_result
