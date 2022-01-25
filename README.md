# DEPair

---

DEPair evaluates ligand-receptor mediated crosstalks among cell populations. 
A interaction score as well as its significance is given for each pair of clusters with DEPair.
DEPair also identifies the variable interactions between different groups.

DEPair is really easy.

```
import numpy as np,scanpy as sc,anndata as ad,pandas as pd
import os,pickle,random
import matplotlib.pyplot as plt
from dfply import *

merged= sc.read_h5ad("./merged.h5ad")
adata_normal      = merged[ merged.obs.status == "normal" ]
adata_degenerated = merged[ merged.obs.status == "abnormal" ]

LRDB=pd.read_table("./CellTalkDB_human_lr_pair.txt")


import depair as dp
from depair.eval import crosstalk

crosstalk_normal=crosstalk(anndata=adata_normal, 
                           LRDB_table=LRDB, 
                           groupby="ann210819",
                           verbose=True, pgb_notebook=True, 
                           n_perm=1000,  n_jobs=48)
                           
crosstalk_degenerated=crosstalk(anndata=adata_degenerated, 
                                LRDB_table=LRDB,
                                groupby="ann210819",
                                verbose=True, pgb_notebook=True, 
                                n_perm=1000,  n_jobs=48)


```


---
## Development plan

- [x] crosstalk intensity evaluation 
- [x] permutation-based crosstalk significance calculation
- [x] permutation speed up
- [x] human and mouse ligand-receptor databases
- [ ] ligand-receptor databases for other species 
- [ ] a tutorial as well as a document online at readthedocs
- [ ] upload to pypi
- [ ] spatial transcriptomics support
- [ ] multi subunit analysis support
- [ ] ligand-receptor database with detailed annotation
- [ ] clustering granularity selection
- [ ] visualization
- [ ] scanpy support
- [ ] supergraph exploration

---

## News
- 2021-04-08  prerelease test version online
- 2021-04-30  Qijin Yin contributed the permutation speed up codes

---

contributors: Sijie Chen, Qijin Yin
