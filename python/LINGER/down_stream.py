import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
Datadir='/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/'
outdir=Datadir + "output/"#output dir
metadata = pd.read_csv((Datadir + "data/metadata.csv"), index_col=0)
metadata['group'] = (metadata['group'] == 'HEROIN').astype(int)
TG_pseudobulk = pd.read_csv((Datadir + "data/TG_pseudobulk.tsv"), index_col=0)
Gene_score = pd.read_csv("/project/zduren/durenlab/palmetto/cham/Heroin/ori_data/GENE_score_uniq.tsv" ,sep = '\t', index_col=0)
K=15
p=0.9
plt.figure(figsize=(10, 6))
def Module_trans(outdir,metadata,TG_pseudobulk,K,GWASfile=None):
    import numpy as np
    from scipy import stats

    
    print('loading GRN......')
    
    trans_reg=pd.read_csv(outdir+'cell_population_trans_regulatory.txt',sep='\t',index_col=0)
    #TG_pseudobulk=TG_pseudobulk_all
    TFset=trans_reg.columns
    TGset=trans_reg.index
    #TG_pseudobulk=TG_pseudobulk/TG_pseudobulk.mean()
    #idx = [s[:3] != "MIR" for s in TG_pseudobulk.index]
    #TG_pseudobulk=TG_pseudobulk.loc[TG_pseudobulk.index[idx]]
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
    print('differential modules......')
    celltype = metadata['celltype'].unique()
    pvalue_all,tvalue_all=diff_Module(Exp_TG,metadata,celltype,S_TG,K)
    nlog_P = -np.log10(pvalue_all)
    sns.heatmap(nlog_P)
    plt.savefig((Datadir + str(K) +"Module_MAXiter.png"),dpi=300)
    plt.clf()

    for Line in Line_list:
        sex_idx = (metadata['male']==1).values
        Line_idx = (metadata['Line']==Line).values
#        celltype_idx = (metadata['celltype']!='46').values
#        celltype_idx = ~(metadata['celltype'].isin([46,47])).values
        pvalue_all,tvalue_all=diff_Module(Exp_TG.loc[:,Line_idx&sex_idx],metadata[Line_idx&sex_idx],celltype,S_TG,K)
     #   pvalue_all,tvalue_all=diff_Module(Exp_TG.loc[:,sex_idx&celltype_idx],metadata[sex_idx],celltype,S_TG,K)
    #    pvalue_all,tvalue_all=diff_Module(Exp_TG.loc[:,sex_idx],metadata[sex_idx],celltype,S_TG,K)

    #    pvalue_all,tvalue_all=diff_Module(Exp_TG,metadata,celltype,S_TG,K)
        nlog_P = -np.log10(pvalue_all+0.001)
    #    nlog_P = -np.log10(pvalue_all)
        sns.heatmap(nlog_P.loc[:,0:46])
        plt.title(Line)
        plt.savefig("./Module_pval_heatmap_M_" + Line + str(K) +".png",dpi=300)
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





#S_TG[(S_TG==1).values].index.values


from scipy.stats import ttest_1samp

def diff_Module_cham(metadata,celltype_list,S_TG,K,DEG_df):
    pvalue_all= np.zeros((K, len(celltype_list)))
    pvalue_all= pd.DataFrame(pvalue_all,index=['M'+str(i+1) for i in range(K)],columns=celltype_list)
    tvalue_all= np.zeros((K, len(celltype_list)))
    tvalue_all= pd.DataFrame(tvalue_all,index=['M'+str(i+1) for i in range(K)],columns=celltype_list)
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
#    for k in range(len(celltype)):
    for celltype in celltype_list:
#        temp=Exp_TG.iloc[:,metadata['celltype'].values==celltype]
 #       aud_idxtemp=metadata[(metadata['celltype'].values==celltype)]['group'].values
        for i in range(K):
            print("Module",i)
            logfc_values = DEG_df.loc[S_TG['Module']==i, celltype]
            stat, p_value = ttest_1samp(logfc_values, 0)
            tvalue_all.loc['M'+str(i+1),celltype]=stat
            pvalue_all.loc['M'+str(i+1),celltype]=p_value
        # Print the results
    return pvalue_all,tvalue_all


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



for i in range(K+1):
     print("Module ",i)
     print((S_TG['Module']==i).sum())




import statsmodels.formula.api as smf
import pandas as pd

def fit_model(data, gene):
    formula = f'{gene} ~ Age + C(Sex) + C(CellType)'  # Adjust formula as needed
    model = smf.ols(formula, data=data).fit()
    return model.pvalues, model.params  # You can adjust what to return based on your needs



def Module_trans(outdir,metadata,TG_pseudobulk,K,p,GWASfile=None):
    import numpy as np
    from scipy import stats

    
    print('loading GRN......')
    
    trans_reg=pd.read_csv(outdir+'cell_population_trans_regulatory.txt',sep='\t',index_col=0)
    trans_reg_cell=pd.read_csv(outdir+ 'cell_type_specific_trans_regulatory_L6_IT_Car3.txt',sep='\t',index_col=0)
    TFset=trans_reg.columns
    TGset=trans_reg.index
    R1 = stats.zscore(trans_reg,1);R1[np.isnan(R1)] = 0.0
    R2 = stats.zscore(trans_reg,0);R2[np.isnan(R2)] = 0.0
    from sklearn.preprocessing import quantile_transform
    Exp=pd.DataFrame(TG_pseudobulk,index=TG_pseudobulk.index,columns=TG_pseudobulk.columns)
    Z=R1+R2
    Z[Z<0]=0  
    use_GWAS = False
    if use_GWAS:
        alpha = 1
        TF_Gene_score = pd.DataFrame((~Z.columns.isin(Gene_score.index)).astype(int), index = Z.columns)
        tmp_Gidx = TF_Gene_score[TF_Gene_score[0]==0].index
        TF_Gene_score[TF_Gene_score[0]==0] = Gene_score.loc[tmp_Gidx]
        TF_Gene_score = (1 - (alpha * (1 + np.log10(TF_Gene_score))))

        TG_Gene_score = pd.DataFrame((~Z.index.isin(Gene_score.index)).astype(int), index = Z.index)
        tmp_Gidx = TG_Gene_score[TG_Gene_score[0]==0].index
        TG_Gene_score[TG_Gene_score[0]==0] = Gene_score.loc[tmp_Gidx]
        TG_Gene_score = (1 - (alpha * (1 + np.log10(TG_Gene_score))))
        #TG_Gene_score = (1 -np.log10(TG_Gene_score))/2

        TF_score_matrix = np.diag(TF_Gene_score.values.flatten())
        TG_score_matrix = np.diag(TG_Gene_score.values.flatten())

        old_Z = Z
        Z = TG_score_matrix @ Z# @ TF_score_matrix


    #figure_dir = Datadir+'figure/celltyep_specific/'
    figure_dir = Datadir+'figure/LINGER_module_compare/'

    celltype = metadata['celltype'].unique()
    print('identify modules......')
    from sklearn.decomposition import NMF


    #fig, axes = plt.subplots(9, 6, figsize=(30, 42))
    #fig.suptitle('Module heatmap')
    l=0
   # for k in range(10,16):
    for k in [7,8,9]:
        K = k
        nmf = NMF(n_components=K, init='random', random_state=0,max_iter=1000)
        W = nmf.fit_transform(Z)
        H = nmf.components_
        #for i in range(50,96,5):
        for i in [70]:
            p = round(i * 0.01, 2)
            [S_TG,W2]=assignLabel(W,p);
            [S_TF,H2]=assignLabel(H.T,p);
            #Exp_mean=stats.zscore(Exp_TG.T).T.groupby(S_TG).mean()
            S_TG=pd.DataFrame(S_TG,index=TGset,columns=['Module'])
            S_TF=pd.DataFrame(S_TF,index=TFset,columns=['Module'])
            Exp_TF=Exp.loc[TFset]
            Exp_TG=Exp.loc[TGset]
            print('differential modules......')
            pvalue_all,tvalue_all=diff_Module(Exp_TG,metadata,celltype,S_TG,K)

            for m in range(1,K+1):
                S_TG[S_TG['Module'] == m].to_csv(Datadir + "/Module_genes/K"+ str(K) +"p" + str(p) +"M"+ str(m)  + ".csv")

            nlog_P = -np.log10(pvalue_all)

            sns.heatmap(nlog_P.T, cmap="viridis")
            plt.title("K"+ str(K) +"p" + str(p))
            #plt.title( "Module_genes/K"+ str(K) +"p" + str(p) +"M"+ str(i))
            plt.savefig((figure_dir + ("K"+ str(K) +"p" + str(p)) + "Module_MAXiter.png"),dpi=300)
            plt.clf()
        #axes[l].set_title("K"+ str(K) +"p" + str(p) +"M"+ str(i))
        #plt.savefig((Datadir+ 'figure/K15/' + str(p) +"GWAS_Module_MAXiter.png"),dpi=300)
        l = l+1



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




for i in range(1,K+1):
    S_TG[S_TG['Module'] == i].to_csv(Datadir"/Module_genes/K"+ str(i) +"p" + str(p) + ".csv")
