import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os,pickle,random
import matplotlib.pyplot as plt
import pickle

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
    from tqdm.autonotebook import trange
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
    
#-----------------------------------------------------------------------------
LRDB=pd.read_csv("./LR_DASkelly.csv")
adata, LRDB = primary_trimming(adata, LRDB, groupby="celltype0627")
crosstalk = get_crosstalk(adata, groupby="celltype0627", LRDB=LRDB)
import warnings;warnings.filterwarnings("ignore");
crosstalk2=get_significance(adata, crosstalk, LRDB, n_perm=1000, groupby="celltype0627")

df = pd.DataFrame()
for key in crosstalk2.keys():
    cf=crosstalk2[key]
    cf=cf[np.logical_and(cf.p_val_adj<0.05, cf.log_avgFC>0.5)]
    df=df.append(cf)
