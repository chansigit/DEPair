import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

def iterator_pgbar(iterator, pgb_notebook, ncols=80 ):
    """
    A util function for automatic tqdm progressbar selection
    """
    def in_notebook():
        """
        Returns ``True`` if the module is running in IPython kernel,
        ``False`` if in IPython shell or other Python shell.
        """
        import sys;return 'ipykernel' in sys.modules;
    
    if in_notebook() and pgb_notebook:
        from tqdm.notebook import tqdm as tqdmnb
        return tqdmnb(iterator)
    else:
        from tqdm import tqdm
        return tqdm(iterator, ncols=80)
    
    

def __grouped_obs_mean(adata, group_key, layer=None, 
                     log_space_averaging= True):
    """
    A util function calculating average expression with in each group.
    adata: an AnnData object
    group_key: a column in adata.obs 
    layer: data layer. if None, adata.X is used
    log_space_averaging: if True, averaging in non-log space, and return the log of the averaged values
    """
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
        if log_space_averaging: # averaging in the non-log-space, return the log of the averaged values
            out[group] = np.ravel(np.log(np.exp(X).mean(axis=0, dtype=np.float64))) # axis=0:= avg across cells
        else:
            out[group] = np.ravel(X.mean(axis=0, dtype=np.float64)) # axis=0:= avg across cells
    return out
    

def __primary_trimming(adata, LRDB, groupby):
    """
    A util function exclude genes that are not present in the LRDB table
    """
    # 1. check if gene symbols are occurred in the dataset
    sel=np.logical_and(LRDB["ligand_symbol"].isin(adata.var_names), 
                       LRDB["receptor_symbol"].isin(adata.var_names)) 
    LRDB=LRDB[sel]

    genesel = np.unique(np.union1d(LRDB["ligand_symbol"], 
                                   LRDB["receptor_symbol"]))
    adata   = adata[:,genesel].copy()
    
    # 2. exclude zero expression pairs 
    avgexpr     = __grouped_obs_mean(adata, group_key=groupby)
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

def __get_expression_fraction(anndata, genes, groupby, cluster):
    """
    given AnnData, a list of genes, and one cluster,
    report the expressing fraction of given genes in that cluster
    """
    sel = anndata.obs[groupby]==cluster  # select cell cluster
    subset = anndata[sel, genes]         # select genes
    expr_cnt = (subset.X >0).sum(axis=0) # binarize genes' expressions
    return pd.DataFrame( (expr_cnt / subset.X.shape[0] ).astype(np.float32),
                          index=genes, columns=["frac"] )


    
def __get_crosstalk(adata, groupby, LRDB, verbose=False, pgb_notebook=False, ncols=80):
    # prepare ligand/receptor expressions
    avgexpr     = __grouped_obs_mean(adata, group_key=groupby)
    lig_symbols = LRDB["ligand_symbol"]
    rec_symbols = LRDB["receptor_symbol"]
    lig_expr    = avgexpr.loc[lig_symbols]
    rec_expr    = avgexpr.loc[rec_symbols]
    pair_symbols= LRDB["pair"]
    
    
    # compute cell-type level crosstalk
    crosstalk_by_celltype=dict()
    if verbose:
        iters = iterator_pgbar( avgexpr.columns, pgb_notebook =pgb_notebook,ncols=ncols)
    else:
        iters =  avgexpr.columns
    for ct1 in iters: # sender cell type
        for ct2 in (avgexpr.columns):  # receiver cell type
            lig_vec = lig_expr[ct1].to_numpy()
            rec_vec = rec_expr[ct2].to_numpy()
            df = pd.DataFrame(data=lig_vec*rec_vec,
                              index=pair_symbols, 
                              columns=["intensity"], 
                              dtype=np.float64 )
            df["sender_celltype"]   = ct1
            df["receiver_celltype"] = ct2
            df["ligand_symbol"]   = lig_symbols.tolist()
            df["receptor_symbol"] = rec_symbols.tolist()
            df["ligand_avgexpr"]  = lig_vec
            df["receptor_avgexpr"]= rec_vec
            
            # compute the expressing fraction within sender/receiver
            df["ligand_frac"]  = __get_expression_fraction(anndata=adata, 
                       groupby=groupby, genes=df["ligand_symbol"],   cluster=ct1).to_numpy()
            df["receptor_frac"]= __get_expression_fraction(anndata=adata, 
                       groupby=groupby, genes=df["receptor_symbol"], cluster=ct2).to_numpy()
            
            
            # placeholder for the significances
            df["p_val"]     = np.nan
            df["p_val_adj"] = np.nan
            
            # update the dataframe
            crosstalk_by_celltype[ (ct1, ct2) ]=df
            
    return(crosstalk_by_celltype)
    

def __grouped_rand_mean(adata, group_key, layer=None,log_space_averaging=True):   
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


def __get_significance(adata, crosstalk, LRDB, n_perm, groupby, 
                     pseudo_expr=1e-5,
                     verbose=False, pgb_notebook=False,ncols=80):
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
    if verbose:
        iters = iterator_pgbar( range(n_perm), pgb_notebook=pgb_notebook,ncols=ncols)
    else:
        iters = range(n_perm)
    
    for r in iters: 
        # generate random ligand/receptor expressions
        randexpr      = __grouped_rand_mean(adata, group_key=groupby)
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
            left_tail=2.0*(1-(nLarger[(i,j)]/(1+n_perm)))
            p_val_vec= np.minimum(
                right_tail,left_tail
            )
            crosstalk[(i,j)]["p_val"] =p_val_vec
            crosstalk[(i,j)]["p_val_adj"] =multipletests(pvals=p_val_vec, 
                                                     alpha=0.05, 
                                                     method="fdr_bh",
                                                     is_sorted=False, 
                                                     returnsorted=False)[1]
            crosstalk[(i,j)]["sender_celltype"]= i
            crosstalk[(i,j)]["receiver_celltype"]= j
            crosstalk[(i,j)]["log2_int_null_ratio"] =np.log2(
                 (pseudo_expr +crosstalk[(i,j)]["intensity"].to_numpy())/ 
                 (pseudo_expr +avg_randintensity[(i,j)])  )
          
                
    return crosstalk


from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from tqdm import tqdm

global __GLOBAL_crosstalk, __GLOBAL_LRDB, __GLOBAL_adata, __GLOBAL_groupby

def __poolmap(function,array, n_jobs=30,use_kwargs=None,use_args=None):
    if use_kwargs is None: use_kwargs = True if isinstance(array[0],dict) else False
    if use_args is None: use_args = True if isinstance(array[0],(list,np.ndarray)) else False
    with ProcessPoolExecutor(max_workers=n_jobs) as pool:
        if use_kwargs: futures = [pool.submit(function, **a) for a in array]
        elif use_args:futures = [pool.submit(function, *a) for a in array]
        else:futures = [pool.submit(function, a) for a in array]      
        _ = list( tqdm(as_completed(futures), total=len(futures), ncols=60))
    out = [ x.result() for x in futures]
    return  out

def __parallel_func(perm_idx, ):
    global __GLOBAL_crosstalk, __GLOBAL_LRDB, __GLOBAL_adata, __GLOBAL_groupby
    
    lig_symbols  = __GLOBAL_LRDB["ligand_symbol"]
    rec_symbols  = __GLOBAL_LRDB["receptor_symbol"]
    randexpr     = __grouped_rand_mean(__GLOBAL_adata, group_key=__GLOBAL_groupby)
    randexpr_lig = randexpr.loc[lig_symbols]
    randexpr_rec = randexpr.loc[rec_symbols]
    num_LRs = __GLOBAL_LRDB.shape[0]
    celltype_list = list(set([x[0] for x in __GLOBAL_crosstalk.keys()]))
    
    num_ct = len(celltype_list)
    avg_randintensity = np.zeros((num_ct,num_ct,num_LRs))
    nLarger = np.zeros((num_ct,num_ct,num_LRs))

    # enumerate cell type pairs
    for idx_i,i in  enumerate(celltype_list): # ligand cell type
        for idx_j,j in enumerate(celltype_list): # receptor cell type
            # calculate pair intensities (real/random)
            intensity     = __GLOBAL_crosstalk[(i,j)]['intensity'].values
            randintensity = np.array(randexpr_lig[i]) * np.array(randexpr_rec[j])

            # calculate the average of random intensities
            avg_randintensity[idx_i,idx_j] = randintensity
            nLarger[idx_i,idx_j] = (intensity>=randintensity).astype(int)

    return avg_randintensity, nLarger


def __get_significance_parallel(adata, crosstalk, LRDB, n_perm, groupby, 
                                pseudo_expr=1e-5, 
                                n_jobs=8):
    lig_symbols = LRDB["ligand_symbol"  ]  # get ligand gene symbols
    rec_symbols = LRDB["receptor_symbol"]  # get receptor gene symbols

    celltype_list = list(set([x[0] for x in crosstalk.keys()])) # get cluster list

    # make data shared among jobs
    global __GLOBAL_crosstalk, __GLOBAL_LRDB, __GLOBAL_adata, __GLOBAL_groupby
    import copy
    __GLOBAL_crosstalk = copy.copy(crosstalk)
    __GLOBAL_LRDB      = copy.copy(LRDB)
    __GLOBAL_adata     = copy.copy(adata)
    __GLOBAL_groupby   = copy.copy(groupby)
    
    # generate random permutations
    pool_res = __poolmap(__parallel_func, range(n_perm), n_jobs = n_jobs)
    pool_res = np.array(pool_res) # [n_perm,2,num_ct,num_ct,num_LRs]
    avg_randintensity = pool_res[:,0].mean(axis=0)
    nLarger  = pool_res[:,1].sum(axis=0)
    
    # extreme event counting
    right_tail_mat = 2.0*     (1+nLarger)/(1+n_perm)
    left_tail_mat  = 2.0* (1- (  nLarger)/(1+n_perm) )
    p_val_vec_mat  = np.minimum(right_tail_mat,left_tail_mat)
    
    # calculate p-values and fold changes
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
            crosstalk[(i,j)]["sender_celltype"]= i
            crosstalk[(i,j)]["receiver_celltype"]= j
            crosstalk[(i,j)]["log2_int_null_ratio"] =np.log2(
                 (pseudo_expr +crosstalk[(i,j)]["intensity"].to_numpy())/ 
                 (pseudo_expr +avg_randintensity[idx_i,idx_j])  )
    
                
    return crosstalk    


def crosstalk(anndata, LRDB_table, groupby, n_perm=1000, n_jobs=1, 
              verbose=False, pgb_notebook=False, pseudo_expr=1e-5):
    """
    crosstalk analysis
    anndata: an AnnData object
    LRDB_table: a pandas.DataFrame contains `ligand_symbol` column and `receptor_symbol` column
    groupby: a column key in adata.obs specifying cell types
    n_perm: number of permutations to calculate crosstalk intensity significance. Default = 1000
    n_jobs: number of parallel processes. Default=1

    """
    assert(n_jobs>0 and isinstance(n_jobs ,int) )
    
    # drop the genes that are not present in the data 
    adata_trim, LRDB = __primary_trimming(anndata, LRDB_table, groupby=groupby)
    
    # get intensity scores of the cellular crosstalks
    crosstalk   = __get_crosstalk(adata_trim, groupby=groupby, LRDB=LRDB, verbose=verbose, pgb_notebook=pgb_notebook)
    
    
    # get the significance of the intensities
    if verbose:
        print("Evaluation of interaction significances: randomly permuting %d times"%n_perm)
    else:
        import warnings; warnings.filterwarnings("ignore");
    
    if n_jobs==1: # serial exec
        if verbose: print("single-core run");
        crosstalk = __get_significance(adata_trim, crosstalk, LRDB, n_perm=n_perm, groupby=groupby,
                                       pseudo_expr=pseudo_expr,
                                       verbose=verbose, pgb_notebook=pgb_notebook )
    else:
        if verbose: print("%d-core run"%n_jobs);
        crosstalk = __get_significance_parallel(adata_trim, crosstalk, LRDB, n_perm=n_perm, groupby=groupby,
                                                pseudo_expr=pseudo_expr,
                                                n_jobs=n_jobs )
    
    return crosstalk
