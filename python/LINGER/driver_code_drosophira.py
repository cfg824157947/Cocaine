import pandas as pd
RE_pseudobulk_all=pd.read_csv('data/RE_pseudobulk.tsv',index_col=0,header=0)
TG_pseudobulk_all=pd.read_csv('data/TG_pseudobulk.tsv',index_col=0,header=0)
K=5
GRN='TF_RE_binding'
adjust_method='bonferroni'
corr_method='pearsonr'



def driver_score(expression,aud_idx,GRN,outdir,adjust_method,corr_method,allcelltype=None):
    print('loading GRN......')
    import numpy as np   
    from statsmodels.stats.multitest import multipletests
    allcelltype=aud_idx['celltype'].unique()
#    allcelltype=['alpha/beta/gamma Kenyon cell','R7/R8','Fat body cells']

    C_result=pd.DataFrame([])
    P_result=pd.DataFrame([])
    Q_result=pd.DataFrame([])
    reg=pd.read_csv(outdir+'cell_population_'+GRN+'.txt',sep='\t',index_col=0)
    #print(reg.shape)
    if reg.shape[1]<4:
        reg = reg.drop_duplicates(subset=['RE', 'TF'])
        reg = reg.pivot(index='RE', columns='TF', values='score').fillna(0)
    reg = reg.fillna(0)
    cols=reg.sum(axis=0).values
    rows=reg.sum(axis=1).values
    E=np.reshape(rows,(rows.shape[0],1))*np.reshape(cols,(1,cols.shape[0]))/rows.sum()
    #print(E.mean().mean()*10**(-4))
    E=E+E.mean().mean()*10**(-4)
    reg=(reg-E)/E
    reg[reg<0]=0
    overlap=list(set(expression.index)&set(reg.index))
    reg=reg.loc[overlap]
    expression=expression.loc[overlap]
    for celltype in allcelltype:
        print('cell type '+ celltype)
        aud_idx1=aud_idx.reset_index()
        celltype_barcode = aud_idx[aud_idx['celltype']==celltype].index
        aud_barcode = aud_idx[aud_idx['group']==1].index
        celltype_aud_barcode = celltype_barcode[celltype_barcode.isin(aud_barcode)]
        celltype_ctl_barcode = celltype_barcode[~celltype_barcode.isin(aud_barcode)]
        Mean1=expression[celltype_aud_barcode].mean(axis=1)+10**(-6)
        Mean0=expression[celltype_ctl_barcode].mean(axis=1)+10**(-6)
        FC=pd.DataFrame(Mean1/Mean0)
        FC=FC.loc[reg.index]
        #print(FC)
        print(np.isnan(FC).sum())
        c,cp=correlation_FC(np.log(FC[0]).fillna(0).values,reg,corr_method)
        #c,cp=correlation_FC(np.log(FC[0]).values,reg,corr_method)
        #idx=pd.DataFrame(range(expression.shape[0]),index=expression.index)
        C_result=pd.concat([C_result,c],axis=1)
        P_result=pd.concat([P_result,cp],axis=1)   
        cp=cp.fillna(1)  
        adjusted_p_values = pd.DataFrame(multipletests(cp[0].values, method=adjust_method)[1],index=c.index)
        #print(adjusted_p_values)
        Q_result=pd.concat([Q_result,adjusted_p_values],axis=1)  
    C_result.columns=allcelltype
    P_result.columns=allcelltype
    Q_result.columns=allcelltype    
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

def correlation_FC(x,y,method):
    import numpy as np
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



def driver_score(expression,aud_idx,GRN,outdir,adjust_method,corr_method,allcelltype=None):
    print('loading GRN......')
    import numpy as np   
    from statsmodels.stats.multitest import multipletests
    allcelltype=aud_idx['celltype'].unique()
#    allcelltype=['alpha/beta/gamma Kenyon cell','R7/R8','Fat body cells']

    C_result=pd.DataFrame([])
    P_result=pd.DataFrame([])
    Q_result=pd.DataFrame([])
    reg=pd.read_csv(outdir+'cell_population_'+GRN+'.txt',sep='\t',index_col=0)
    #print(reg.shape)
    if reg.shape[1]<4:
        reg = reg.drop_duplicates(subset=['RE', 'TF'])
        reg = reg.pivot(index='RE', columns='TF', values='score').fillna(0)
    reg = reg.fillna(0)
    cols=reg.sum(axis=0).values
    rows=reg.sum(axis=1).values
    E=np.reshape(rows,(rows.shape[0],1))*np.reshape(cols,(1,cols.shape[0]))/rows.sum()
    #print(E.mean().mean()*10**(-4))
    E=E+E.mean().mean()*10**(-4)
    reg=(reg-E)/E
    reg[reg<0]=0
    overlap=list(set(expression.index)&set(reg.index))
    reg=reg.loc[overlap]
    expression=expression.loc[overlap]
    #cols=expression.sum(axis=0).values
    #rows=expression.sum(axis=1).values
    #E=np.reshape(rows,(rows.shape[0],1))*np.reshape(cols,(1,cols.shape[0]))/rows.sum()
    #E=E+E.mean().mean()*10**(-1)
    #expression=expression/E
    for i in range(len(allcelltype)):
        print('cell type '+ allcelltype[i])
        aud_idx1=aud_idx.reset_index()
        exp_temp=expression.iloc[:,aud_idx1[aud_idx1['celltype']==allcelltype[i]].index]
        aud_idx1=aud_idx1[aud_idx1['celltype']==allcelltype[i]]
        Mean1=exp_temp[aud_idx1[aud_idx1['group']==1]['index']].mean(axis=1)+10**(-6)
        Mean0=exp_temp[aud_idx1[aud_idx1['group']==0]['index']].mean(axis=1)+10**(-6)
        FC=pd.DataFrame(Mean1/Mean0)
        FC=FC.loc[reg.index]
        #print(FC)
        print(np.isnan(FC).sum())
        c,cp=correlation_FC(np.log(FC[0]).fillna(0).values,reg,corr_method)
        #c,cp=correlation_FC(np.log(FC[0]).values,reg,corr_method)
        #idx=pd.DataFrame(range(expression.shape[0]),index=expression.index)
        C_result=pd.concat([C_result,c],axis=1)
        P_result=pd.concat([P_result,cp],axis=1)   
        cp=cp.fillna(1)  
        adjusted_p_values = pd.DataFrame(multipletests(cp[0].values, method=adjust_method)[1],index=c.index)
        #print(adjusted_p_values)
        Q_result=pd.concat([Q_result,adjusted_p_values],axis=1)  
    C_result.columns=allcelltype
    P_result.columns=allcelltype
    Q_result.columns=allcelltype    
    return C_result,P_result,Q_result




def driver_score(expression,aud_idx,GRN,outdir,adjust_method,corr_method,allcelltype=None):
    print('loading GRN......')
    from statsmodels.stats.multitest import multipletests
    allcelltype=aud_idx['celltype'].unique()
#    allcelltype=['alpha/beta/gamma Kenyon cell','R7/R8','Fat body cells']

    C_result=pd.DataFrame([])
    P_result=pd.DataFrame([])
    Q_result=pd.DataFrame([])
    reg=pd.read_csv(outdir+'cell_population_'+GRN+'.txt',sep='\t',index_col=0)
    #print(reg.shape)
    if reg.shape[1]<4:
        reg = reg.drop_duplicates(subset=['RE', 'TF'])
        reg = reg.pivot(index='RE', columns='TF', values='score').fillna(0)
    reg = reg.fillna(0)
    cols=reg.sum(axis=0).values
    rows=reg.sum(axis=1).values
    E=np.reshape(rows,(rows.shape[0],1))*np.reshape(cols,(1,cols.shape[0]))/rows.sum()
    #print(E.mean().mean()*10**(-4))
    E=E+E.mean().mean()*10**(-4)
    reg=(reg-E)/E
    reg[reg<0]=0
    overlap=list(set(expression.index)&set(reg.index))
    reg=reg.loc[overlap]
    expression=expression.loc[overlap]
    for celltype in allcelltype:
        print('cell type '+ celltype)
        aud_idx1=aud_idx.reset_index()
        celltype_barcode = aud_idx[aud_idx['celltype']==celltype].index
        aud_barcode = aud_idx[aud_idx['group']==1].index
        celltype_aud_barcode = celltype_barcode[celltype_barcode.isin(aud_barcode)]
        celltype_ctl_barcode = celltype_barcode[~celltype_barcode.isin(aud_barcode)]
        Mean1=expression[celltype_aud_barcode].mean(axis=1)+10**(-6)
        Mean0=expression[celltype_ctl_barcode].mean(axis=1)+10**(-6)
        FC=pd.DataFrame(Mean1/Mean0)
        FC=FC.loc[reg.index]
        #print(FC)
        print(np.isnan(FC).sum())
        c,cp=correlation_FC(np.log(FC[0]).fillna(0).values,reg,corr_method)
        #c,cp=correlation_FC(np.log(FC[0]).values,reg,corr_method)
        #idx=pd.DataFrame(range(expression.shape[0]),index=expression.index)
        C_result=pd.concat([C_result,c],axis=1)
        P_result=pd.concat([P_result,cp],axis=1)   
        cp=cp.fillna(1)  
        adjusted_p_values = pd.DataFrame(multipletests(cp[0].values, method=adjust_method)[1],index=c.index)
        #print(adjusted_p_values)
        Q_result=pd.concat([Q_result,adjusted_p_values],axis=1)  
    C_result.columns=allcelltype
    P_result.columns=allcelltype
    Q_result.columns=allcelltype    
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


GRN='trans_regulatory'
for i in range(int(np.round(len(cell_types)/5))):
    cells = cell_types[i*5:i*5+5]
    C_result_RNA,P_result_RNA,Q_result_RNA=driver_score(TG_pseudobulk_all,metadata,GRN,outdir,adjust_method,corr_method,cells)
    C_result_TG_r,Q_result_TG_r=driver_result(C_result_RNA,Q_result_RNA,K)
    Q_result_TG_r=Q_result_TG_r.clip(lower=10**-100)
    C_result_TG_r.to_csv('C_result_TG_r'+str(i)+'.txt',sep='\t')
    Q_result_TG_r.to_csv('Q_result_TG_r'+str(i)+'.txt',sep='\t')
    print(cells)




# Read the data
dataP = pd.read_csv('Q_result_TG_r.txt', sep='\t', index_col=0)
dataT = pd.read_csv('C_result_TG_r.txt', sep='\t', index_col=0)

# Transform dataP
dataP = -np.log10(dataP)
dataP = dataP.clip(upper=100)  # Cap values at 100

# Add TF as a column for melting
dataP['TF'] = dataP.index
dataT['TF'] = dataT.index

# Melt data for long format
longdiff0 = pd.melt(dataP, id_vars=['TF'], var_name='celltype', value_name='P')
longdiff1 = pd.melt(dataT, id_vars=['TF'], var_name='celltype', value_name='PCC')

# Merge the data
longdiff1['P'] = longdiff0['P']

# Define plot parameters
cutoff = 50
maxp = np.ceil(longdiff1['P'].max())
limits0 = (cutoff, maxp)
range0 = (1, 4)
numbreak = 5
d = np.ceil((maxp - cutoff) / (numbreak - 1))
breaks0 = np.linspace(cutoff, cutoff + d * (numbreak - 1), numbreak)

# Plot
plt.figure(figsize=(10, 8))
scatter = plt.scatter(
    x=longdiff1['celltype'],
    y=longdiff1['TF'],
    s=longdiff1['P'].map(lambda x: range0[0] + (x - cutoff) / (limits0[1] - limits0[0]) * (range0[1] - range0[0])),
    c=longdiff1['PCC'],
    cmap='RdBu_r',
    edgecolors='k',
    alpha=1
)

# Colorbar and size legend
plt.colorbar(scatter, label='PCC')
plt.xticks(rotation=90, fontsize=9, ha='right')
plt.yticks(fontsize=10)
plt.xlabel("Cell Type")
plt.ylabel("TF")

# Size legend
handles = [plt.scatter([], [], s=size, edgecolors='k', color='gray') for size in range0]
labels = [f"{b:.1f}" for b in breaks0]
plt.legend(handles, labels, title="P", labelspacing=1.2, loc='right', bbox_to_anchor=(1.2, 0.5))




r('''
library(ggplot2)
library(grid)

for(iter in c(0:8)){
P_file = paste0('Q_result_TG_r',i,'.txt')
T_file = paste0('C_result_TG_r',i,'.txt')
dataP=read.table(P_file,sep='\t',row.names=1,header=TRUE)
dataT=read.table(T_file,sep='\t',row.names=1,header=TRUE)
sort_TF=rownames(dataT)
library(tidyr)
dataP=-log10(dataP)
print(paste0('maxinum of -log10P:',max(dataP)))
maxP=100
dataP[dataP>100]=100
dataP1=dataP
dataP1$TF=rownames(dataP)
longdiff0 <- gather(dataP1, sample, value,-TF)
longdiff0_s <- longdiff0[order(longdiff0$TF, longdiff0$sample), ]
dataT1=dataT
dataT1$TF=rownames(dataT)
longdiff1=gather(dataT1, sample, value,-TF)
longdiff1=longdiff1[order(longdiff1$TF, longdiff1$sample), ]
colnames(longdiff1)=c('TF','celltype','PCC')
longdiff1$P=longdiff0_s$value
longdiff1$TF=factor(longdiff1$TF,levels=(sort_TF))
library(egg)
cutoff=50
maxp=ceiling(max(dataP))
print(maxp)
limits0=c(cutoff,maxp)
print(limits0)
range0 = c(1,4)
numbreak=5
d=ceiling((maxp-cutoff)/(numbreak-1))
print(d)
breaks0= seq(from = cutoff, to = cutoff+d*(numbreak-1),  length.out=numbreak)

print(range0)
print(breaks0)
p=ggplot(longdiff1,aes(x = celltype, y = TF))+
geom_point(aes(size = P, fill = PCC), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = limits0, range = range0, breaks = breaks0) + 
  labs( x= "cell type", y = "TF", fill = "")  + theme_article()+
  theme(legend.key=element_blank(), 
  axis.text.x = element_text( size = 9, face = "bold", angle = 0, vjust = 0.3, hjust = 1), 
  legend.position = "right") + 
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",
                     high="red", space ="Lab" )
fileR=paste0("plot",i,".png")
ggsave(filename = "plot.png", plot = p, width = 6, height = 4, dpi = 300)
}
''')
