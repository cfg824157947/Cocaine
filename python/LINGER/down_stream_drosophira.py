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


