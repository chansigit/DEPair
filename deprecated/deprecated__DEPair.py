import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os,pickle,random
import matplotlib.pyplot as plt
import pickle

####  multi-processing
from tqdm import tqdm
#from tqdm.notebook import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
def poolmap(function,array, n_jobs=30,use_kwargs=None,use_args=None):
    if use_kwargs is None: use_kwargs = True if isinstance(array[0],dict) else False
    if use_args is None: use_args = True if isinstance(array[0],(list,np.ndarray)) else False
    with ProcessPoolExecutor(max_workers=n_jobs) as pool:
        if use_kwargs: futures = [pool.submit(function, **a) for a in array]
        elif use_args:futures = [pool.submit(function, *a) for a in array]
        else:futures = [pool.submit(function, a) for a in array]      
        _ = list( tqdm(as_completed(futures), total=len(futures)))
    out = [ x.result() for x in futures]
    return  out
####  multi-processing


def grouped_obs_mean(adata, group_key, layer=None, 
                     log_space_averaging= True):
    # define an operator for expression matrix retieval
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X   

    # split metadata dataframe into sets of rows by group_key
    grouped = adata.obs.groupby(group_key)
    
    # build new empty expression matrix
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64), # nGene-by-nCelltype zero matrix
        columns = list(grouped.groups.keys()), # colnames := list of cell types
        index   = adata.var_names #rowname:= gene symbols
    )
    
    # compute average expression in the celltype-wise manner (by celltype column)
    for group, idx in grouped.indices.items(): # group:=celltype idx:=barcodes
        X = getX(adata[idx]) # cell-by-gene matrix
        if log_space_averaging: # averaging in the log-space, return the log of the averaged values
            out[group] = np.ravel(np.log(np.exp(X).mean(axis=0, dtype=np.float64))) # axis=0:= avg across cells
        else:
            out[group] = np.ravel(X.mean(axis=0, dtype=np.float64)) # axis=0:= avg across cells
    return out
    


def primary_trimming(adata, LRDB, groupby):
    # 1. check if gene symbols are occurred in the dataset
    sel=np.logical_and(LRDB["ligand_symbol"].isin(adata.var_names), 
                       LRDB["receptor_symbol"].isin(adata.var_names)) 
    LRDB=LRDB[sel]

    genesel = np.unique(np.union1d(LRDB["ligand_symbol"], 
                                   LRDB["receptor_symbol"]))
    adata   = adata[:,genesel].copy()
    
    # 2. exclude zero expression pairs 
    avgexpr     = grouped_obs_mean(adata, group_key=groupby)
    lig_symbols = LRDB["ligand_symbol"]
    rec_symbols = LRDB["receptor_symbol"]
    lig_expr    = avgexpr.loc[lig_symbols]
    rec_expr    = avgexpr.loc[rec_symbols]
    sel = np.logical_not(np.logical_or(
              (rec_expr.T==0).all().to_numpy(),
              (lig_expr.T==0).all().to_numpy() ))
    LRDB = LRDB[sel]
    
    genesel = np.unique(np.union1d(LRDB["ligand_symbol"], 
                                   LRDB["receptor_symbol"]))
    adata   = adata[:,genesel].copy()
    
    return(adata, LRDB)


def get_crosstalk(adata, groupby, LRDB):
    # prepare ligand/receptor expressions
    avgexpr     = grouped_obs_mean(adata, group_key=groupby)
    lig_symbols = LRDB["ligand_symbol"]
    rec_symbols = LRDB["receptor_symbol"]
    lig_expr    = avgexpr.loc[lig_symbols]
    rec_expr    = avgexpr.loc[rec_symbols]
    pair_symbols= LRDB["pair"]
    
    # compute cell-type level crosstalk
    crosstalk_by_celltype=dict()
    for ct1 in (avgexpr.columns):
        for ct2 in (avgexpr.columns):
            lig_vec = lig_expr[ct1].to_numpy()
            rec_vec = rec_expr[ct2].to_numpy()
            df = pd.DataFrame(data=lig_vec*rec_vec,
                              index=pair_symbols, 
                              columns=["intensity"], 
                              dtype=np.float64 )
            df["ligand_symbol"]   = lig_symbols.tolist()
            df["receptor_symbol"] = rec_symbols.tolist()
            df["ligand_avgexpr"]  = lig_vec
            df["receptor_avgexpr"]= rec_vec
            df["p_val"]     = np.nan
            df["p_val_adj"] = np.nan
            crosstalk_by_celltype[ (ct1, ct2) ]=df
    return(crosstalk_by_celltype)
    

def grouped_rand_mean(adata, group_key, layer=None,log_space_averaging=True):   
    # define an operator for expression matrix retieval
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    
    # shuffle the identities 
    import random
    adata.obs["random"]=random.sample(adata.obs[group_key].tolist(), len(adata.obs))
    
    # split metadata dataframe into sets of rows by group_key
    grouped = adata.obs.groupby("random")
    
    # build new empty expression matrix
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64), # nGene-by-nCelltype zero matrix
        columns=list(grouped.groups.keys()), # colnames := list of cell types
        index=adata.var_names #rowname:= gene symbols
    )

    # compute average expression in the celltype-wise manner (by celltype column)
    for group, idx in grouped.indices.items(): # group:=celltype idx:=barcodes
        X = getX(adata[idx]) # cell-by-gene matrix
        if log_space_averaging: # averaging in the log-space, return the log of the averaged values
            out[group] = np.ravel(np.log(np.exp(X).mean(axis=0, dtype=np.float64))) # axis=0:= avg across cells
        else:
            out[group] = np.ravel(X.mean(axis=0, dtype=np.float64)) # axis=0:= avg across cells
    return out


def get_significance(adata, crosstalk, LRDB, n_perm, groupby, 
                     use_jupyter=True, pseudo_expr=1e-5):
    lig_symbols = LRDB["ligand_symbol"]
    rec_symbols = LRDB["receptor_symbol"]
    nLarger = dict()
    #nLower  = dict()
    avg_randintensity =dict()
    celltype_list = list(set([x[0] for x in crosstalk.keys()]))
    #celltype_list = list(randexpr.columns)
    
    # generate random permutations
    from tqdm import tqdm
    from tqdm import trange
    for r in trange(n_perm): 
        # generate random ligand/receptor expressions
        randexpr      = grouped_rand_mean(adata, group_key=groupby)
        randexpr_lig  = randexpr.loc[lig_symbols]
        randexpr_rec  = randexpr.loc[rec_symbols]
        
        
        # enumerate cell type pairs
        for i in celltype_list: # ligand cell type
            for j in celltype_list: # receptor cell type
                # calculate pair intensities (real/random)
                intensity     = crosstalk[(i,j)]['intensity']
                randintensity = np.array(randexpr_lig[i]) * np.array(randexpr_rec[j])
                
                # calculate the average of random intensities
                if (i,j) not in avg_randintensity:
                    avg_randintensity[(i,j)] = randintensity
                else:
                    avg_randintensity[(i,j)] = avg_randintensity[(i,j)]+ randintensity
                
                # count extreme events
                if (i,j) not in nLarger:
                    #import pdb;pdb.set_trace();
                    nLarger[(i,j)] = (intensity>=randintensity).astype(int)
                else:
                    nLarger[(i,j)] = nLarger[(i,j)] + (intensity>=randintensity).astype(int)
                #if (i,j) not in nLower:
                #    nLower[(i,j)]  = (intensity<randintensity).astype(int)
                #else:
                #    nLower[(i,j)]  = nLower[(i,j)]  + (intensity<randintensity).astype(int)
    
    # calculate p-values and fold changes
    from statsmodels.stats.multitest import multipletests
    pseudo_expr = 0 if pseudo_expr is None else pseudo_expr
    for i in (celltype_list):
        for j in (celltype_list):
            avg_randintensity[(i,j)] = avg_randintensity[(i,j)]/n_perm
            crosstalk[(i,j)]["avg_null"] = avg_randintensity[(i,j)]
            right_tail=2.0*(1+nLarger[(i,j)])/(1+n_perm)
            #left_tail=2.0*(1+nLower[(i,j)])/(1+n_perm)
            left_tail=2.0*(1-(nLarger[(i,j)]/(1+n_perm)))
            #import pdb;pdb.set_trace();
            p_val_vec= np.minimum(
                right_tail,left_tail
            )
            crosstalk[(i,j)]["p_val"] =p_val_vec
            crosstalk[(i,j)]["p_val_adj"] =multipletests(pvals=p_val_vec, 
                                                     alpha=0.05, 
                                                     method="fdr_bh",
                                                     is_sorted=False, 
                                                     returnsorted=False)[1]
            crosstalk[(i,j)]["cell_ligand"]= i
            crosstalk[(i,j)]["cell_receptor"]= j
            crosstalk[(i,j)]["log_avgFC"] =np.log2(
                 (pseudo_expr +crosstalk[(i,j)]["intensity"].to_numpy())/ 
                 (pseudo_expr +avg_randintensity[(i,j)])  )

            # filter analysis
            #sel = crosstalk[(i,j)].log_avgFC>0
            #crosstalk[(i,j)]=crosstalk[(i,j)][sel]           
                
    return crosstalk            
  

def parallel_func(perm_idx,):
    global GLOBAL_crosstalk,GLOBAL_crosstalk,GLOBAL_LRDB,GLOBAL_groupby
    
    lig_symbols = GLOBAL_LRDB["ligand_symbol"]
    rec_symbols = GLOBAL_LRDB["receptor_symbol"]
    randexpr      = grouped_rand_mean(GLOBAL_adata, group_key=GLOBAL_groupby)
    randexpr_lig  = randexpr.loc[lig_symbols]
    randexpr_rec  = randexpr.loc[rec_symbols]
    num_LRs = GLOBAL_LRDB.shape[0]
    celltype_list = list(set([x[0] for x in GLOBAL_crosstalk.keys()]))
    
    num_ct = len(celltype_list)
    avg_randintensity = np.zeros((num_ct,num_ct,num_LRs))
    nLarger = np.zeros((num_ct,num_ct,num_LRs))
    # enumerate cell type pairs
    for idx_i,i in  enumerate(celltype_list): # ligand cell type
        for idx_j,j in enumerate(celltype_list): # receptor cell type
            # calculate pair intensities (real/random)
            intensity     = GLOBAL_crosstalk[(i,j)]['intensity'].values
            randintensity = np.array(randexpr_lig[i]) * np.array(randexpr_rec[j])

            # calculate the average of random intensities
            avg_randintensity[idx_i,idx_j] = randintensity
            nLarger[idx_i,idx_j] = (intensity>=randintensity).astype(int)
    return avg_randintensity,nLarger

global GLOBAL_crosstalk,GLOBAL_LRDB,GLOBAL_adata,GLOBAL_groupby

def get_significance_parallel(adata, crosstalk, LRDB, n_perm, groupby, 
                     use_jupyter=True, pseudo_expr=1e-5,n_jobs=30):
    lig_symbols = LRDB["ligand_symbol"]
    rec_symbols = LRDB["receptor_symbol"]

    celltype_list = list(set([x[0] for x in crosstalk.keys()]))

    global GLOBAL_crosstalk,GLOBAL_LRDB,GLOBAL_adata,GLOBAL_groupby
    GLOBAL_crosstalk = crosstalk
    GLOBAL_LRDB = LRDB
    GLOBAL_adata = adata
    GLOBAL_groupby = groupby
    
    
    pool_res = poolmap(parallel_func,range(n_perm),n_jobs = n_jobs)
    pool_res = np.array(pool_res) # [n_perm,2,num_ct,num_ct,num_LRs]
    avg_randintensity = pool_res[:,0].mean(axis=0)
    nLarger = pool_res[:,1].sum(axis=0)
    
    #right_tail_mat = 2.0*(1+nLarger)/(1+n_perm) 
    #left_tail_mat = 2.0*(1-nLarger)/(1+n_perm)    # this is a bug, causing negative p-values, corrected as following two lines
    right_tail_mat = 2.0*     (1+nLarger)/(1+n_perm)
    left_tail_mat  = 2.0* (1- (  nLarger)/(1+n_perm) )
    p_val_vec_mat = np.minimum(right_tail_mat,left_tail_mat)
    
    from statsmodels.stats.multitest import multipletests
    pseudo_expr = 0 if pseudo_expr is None else pseudo_expr
    for idx_i,i in enumerate(celltype_list):
        for idx_j,j in enumerate(celltype_list):
            p_val_vec =p_val_vec_mat[idx_i,idx_j]
            
            crosstalk[(i,j)]["avg_null"] = avg_randintensity[idx_i,idx_j]
            crosstalk[(i,j)]["p_val"] =p_val_vec
            crosstalk[(i,j)]["p_val_adj"] =multipletests(pvals=p_val_vec, 
                                                     alpha=0.05, 
                                                     method="fdr_bh",
                                                     is_sorted=False, 
                                                     returnsorted=False)[1]
            crosstalk[(i,j)]["cell_ligand"]= i
            crosstalk[(i,j)]["cell_receptor"]= j
            crosstalk[(i,j)]["log_avgFC"] =np.log2(
                 (pseudo_expr +crosstalk[(i,j)]["intensity"].to_numpy())/ 
                 (pseudo_expr +avg_randintensity[idx_i,idx_j])  )
    
                
    return crosstalk    

#-----------------------------------------------------------------------------

def get_expression_fraction(anndata, gene, groupby, group):
    sel=anndata.obs[groupby]==group
    subset = anndata[sel,:]
    tot=subset.n_obs
    expressed=subset[subset[:, gene].X > 0, :].n_obs
    return float(expressed/tot)

def get_expressing_fraction(adata, genes, groups, groupby):
    gene_ids = list(set(genes))# adata.var.index.values
    clusters = adata.obs[groupby].unique().tolist()
    obs = adata[:,gene_ids].X.toarray()
    obs = pd.DataFrame(obs,columns=gene_ids,index=adata.obs[groupby])
    #average_obs = obs.groupby(level=0).mean()
    obs_bool = obs.astype(bool)
    fraction_obs = obs_bool.groupby(level=0).sum()/obs_bool.groupby(level=0).count()
    #return fraction_obs.T
    #i = pd.Series(group)
    return fraction_obs.T.lookup(genes, groups)

def crosstalk_analysis(anndata, LRDB_table, groupby, n_perm=1000):
    adata, LRDB = primary_trimming(anndata, LRDB_table, groupby=groupby)
    crosstalk = get_crosstalk(adata, groupby=groupby, LRDB=LRDB)
    import warnings;warnings.filterwarnings("ignore");
    crosstalk = get_significance(adata, crosstalk, LRDB, n_perm=n_perm, groupby=groupby)
    
    df_crosstalk = pd.DataFrame()
    for key in crosstalk.keys():
        cf=crosstalk[key]
        df_crosstalk=df_crosstalk.append(cf)
    df_crosstalk.index+="_"+df_crosstalk["cell_ligand"]+"_"
    df_crosstalk.index+=    df_crosstalk["cell_receptor"]
    
    # add gene expression fraction
    ##def add_lig_fraction(row):
    ##    return get_expression_fraction(anndata, row["ligand_symbol"], groupby=groupby, group=row["cell_ligand"])
    ##def add_rec_fraction(row):  
    ##    return get_expression_fraction(anndata, row["receptor_symbol"],groupby=groupby,group=row["cell_receptor"])
    ##df_crosstalk['ligand_pct']  =df_crosstalk.apply (lambda row: add_lig_fraction(row), axis=1)
    ##df_crosstalk['receptor_pct']=df_crosstalk.apply (lambda row: add_rec_fraction(row), axis=1)
    df_crosstalk['ligand_pct']  = get_expressing_fraction(adata=anndata, 
                                                          genes=df_crosstalk["ligand_symbol"],
                                                          groups=df_crosstalk["cell_ligand"],
                                                          groupby=groupby)
    df_crosstalk['receptor_pct']= get_expressing_fraction(adata=anndata, 
                                                          genes=df_crosstalk["receptor_symbol"],
                                                          groups=df_crosstalk["cell_receptor"],
                                                          groupby=groupby)
    return(df_crosstalk)

#---------------------------------------------------------------------------------


def diff_pair(crosstalk1, crosstalk2, 
              ident1, ident2):
    # get common, group1-only, group2-only pairs
    common_idx =  crosstalk1.index.intersection(crosstalk2.index)
    g1_idx = crosstalk1.index.difference(crosstalk2.index)
    g2_idx = crosstalk2.index.difference(crosstalk1.index)
    
    # chosen columns
    selected_col=["intensity","ligand_pct","receptor_pct"]
    
    # assembly group comparison DataFrames for common pairs
    df_group1= crosstalk1.loc[common_idx][selected_col]
    df_group1.rename(columns =  dict(zip(df_group1.columns, df_group1.columns+"_"+ident1)), inplace=True  )
    df_group2= crosstalk2.loc[common_idx][selected_col]
    df_group2.rename(columns =  dict(zip(df_group2.columns, df_group2.columns+"_"+ident2)), inplace=True  )
    df_common= df_group1.join(df_group2)
    
    # assembly group comparison DataFrame for group1-only pairs
    import copy
    df_group1 =crosstalk1.loc[g1_idx][selected_col]
    df_group2 =copy.copy(df_group1) 
    for col in df_group2.columns: 
        df_group2[col].values[:] = 0
    df_group1.rename(columns =  dict(zip(df_group1.columns, df_group1.columns+"_"+ident1)), inplace=True  )
    df_group2.rename(columns =  dict(zip(df_group2.columns, df_group2.columns+"_"+ident2)), inplace=True  )
    df_group1only= df_group1.join(df_group2)
    
    # assembly group comparison DataFrame for group2-only pairs
    df_group2 =crosstalk2.loc[g2_idx][selected_col]
    df_group1 =copy.copy(df_group2)
    for col in df_group1.columns: 
        df_group1[col].values[:] = 0
    df_group1.rename(columns =  dict(zip(df_group1.columns, df_group1.columns+"_"+ident1)), inplace=True  )
    df_group2.rename(columns =  dict(zip(df_group2.columns, df_group2.columns+"_"+ident2)), inplace=True  )
    df_group2only= df_group1.join(df_group2)
    
    df= df_common.append(df_group1only).append(df_group2only)
    
    intensity1=df["intensity_"+ident1]
    intensity2=df["intensity_"+ident2]
    df["intensity_log2FC"]=np.log2(intensity2/intensity1)
    df= df.reindex(sorted(df.columns), axis=1)
    return df

#--------------------------------------------------------------------------------
def depair_plot(DEpair, ident1, ident2, figsize=[7,6],
                xcut=0.5, ycut=0.5, cutwidth=0.5, ptsize=10, vmin=-1, vmax=1,
                highlight=True, highlight_fc_upper=0.5, highlight_fc_lower=-0.5
               ):
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    cmap = sns.color_palette("coolwarm", as_cmap=True,)
    f, ax = plt.subplots(figsize=figsize)
    
    # axis names
    xaxis,yaxis="intensity_"+ident1, "intensity_"+ident2
    ax.set(xlabel="Intensities "+ident1, ylabel="Intensities "+ident2)
    
    # get coordinates
    xpos=DEpair[xaxis]
    ypos=DEpair[yaxis]
    points = ax.scatter(x=xpos, 
                        y=ypos, 
                        c=DEpair["intensity_log2FC"], s=ptsize, cmap=cmap,
                        vmin=vmin,vmax=vmax
                        )
    # highlight points
    if highlight:
        label_sel = np.logical_or(DEpair["intensity_log2FC"]>highlight_fc_upper,
                                  DEpair["intensity_log2FC"]<highlight_fc_lower)
        label_sel = np.logical_and(label_sel, DEpair[xaxis]>xcut)
        label_sel = np.logical_and(label_sel, DEpair[yaxis]>ycut)
        
        label_x = DEpair[label_sel][xaxis].tolist()
        label_y = DEpair[label_sel][yaxis].tolist()
        label_text = DEpair[label_sel].index.tolist()
        
        texts = []
        for x, y, s in zip(label_x, label_y, label_text):
            texts.append(ax.text(x, y, s))
        from adjustText import adjust_text
        adjust_text(texts, force_points=0.2, force_text=0.2,
                    expand_points=(1, 1), expand_text=(1, 1),
                    arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    
    # bound of axis
    bound= DEpair[xaxis].max()
    bound= max(bound, DEpair[yaxis].max())
    bound = int(bound+1)
    
    # set ticks
    ax.xaxis.set_ticks(np.linspace(start=0,stop=bound, num=int(bound+1)))
    ax.yaxis.set_ticks(np.linspace(start=0,stop=bound, num=int(bound+1)))
    from matplotlib.ticker import AutoMinorLocator
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    # draw auxilary lines
    plt.plot([0, bound], [0, bound], 'r--',linewidth=cutwidth)
    plt.plot([xcut, ycut], [0, bound], 'r--',linewidth=cutwidth)
    plt.plot([0,bound], [xcut, ycut], 'r--',linewidth=cutwidth)
    f.colorbar(points)
