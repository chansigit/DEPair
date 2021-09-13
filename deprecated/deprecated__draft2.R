library(Seurat)
library(dplyr)
library(tictoc)
library(ggpubr)

TrimDB.primary<-function(seu, LRDB){
    message("Triming LR-pair DB")
    # trim ligand-receptor pair database 
    ligand.indata   <- unique(LRDB$ligand_symbol[LRDB$ligand_symbol %in%rownames(seu)])
    receptor.indata <- unique(LRDB$receptor_symbol[LRDB$receptor_symbol %in%rownames(seu)])
    lrpair.indata   <-(LRDB$ligand_symbol %in% ligand.indata)&(LRDB$receptor_symbol %in% receptor.indata)
    LRDB            <- LRDB[lrpair.indata, ]
    return(LRDB)
}

TrimSeuratFeatures.primary <- function(seu, LRDB){
    message("Triming Seurat object")
    # trim seurat object
    DefaultAssay(seu) <- "RNA"
    features.sel<-unique(c(LRDB$ligand_symbol, LRDB$receptor_symbol))
    seu<-subset(seu, features = features.sel)
    
    for (assay_name in names(seu@assays)){
        if (assay_name!="RNA"){ 
            seu[[assay_name]]<-NULL
        }
    }
    
    return(seu)
}

TrimDB.secondary <- function(crosstalk.df, LRDB){
    lrpair.detected<-unique(crosstalk.df$pair)
    LRDB<-LRDB[LRDB$pair %in% lrpair.detected, ]
    return(LRDB)
}
TrimSeuratFeatures.secondary<- function(seu, crosstalk.df){
    features.sel<-unique(c(crosstalk.df$ligand, crosstalk.df$receptor))
    seu<-subset(seu, features = features.sel)
    return(seu)
}

GetCrosstalks <- function(seu, LRDB, ident.name, expr.cutoff){
    Idents(seu)<-ident.name
    avgexpr <- AverageExpression(object=seu, assays="RNA",slot="data",
                                 verbose = F, return.seurat = FALSE)$RNA
    df <-data.frame(celltype.1=character(),
                    celltype.2=character(),
                    pair=character(), ligand=character(), receptor=character(),
                    intensity=numeric()
                    )
  
    for (ct1 in levels(seu)){
        for (ct2 in levels(seu)){
            for (pairID in 1:nrow(LRDB)){
                ligand  <-LRDB[pairID, "ligand_symbol"  ]
                receptor<-LRDB[pairID, "receptor_symbol"]
                pair    <-LRDB[pairID, "pair"]
                if (avgexpr[ligand, ct1]>expr.cutoff && 
                    avgexpr[receptor, ct2]>expr.cutoff){
                    intensity <- avgexpr[ligand, ct1]*  avgexpr[receptor, ct2]
                    onerow    <- list(celltype.1=ct1,
                                      celltype.2=ct2, 
                                      pair=pair,
                                      ligand=ligand, receptor=receptor,
                                      intensity=intensity)
                    df<- rbind(df, onerow)
                    #print(paste(ct1,"-> ⊃-",ct2, '\tvia\t',
                    #      ligand," & ",receptor, "\t", intensity))
                }
            }
        }
    }
    return (df)
}

GetRandExpression<-function(rand.batch, seu, gene,celltype){
    rand.ident.name<-paste("random", rand.batch, sep="_")
    Idents(seu)<-rand.ident.name
    avgrandexpr <- AverageExpression(object=subset(seu,idents=celltype),
                                 assays="RNA",slot="data",features = gene,
                                 verbose = F, return.seurat = FALSE)$RNA[gene,]
    return(avgrandexpr)
}





tic();AverageExpression(subset(immune.combined, idents="2.CD4+CTLA4+ T"),assays="RNA", features = "Ltf",verbose = F)$RNA[,"2.CD4+CTLA4+ T"];toc()
tic();AverageExpression(immune.combined,assays="RNA", features = "Ltf",verbose = F)$RNA[,"2.CD4+CTLA4+ T"];toc()



ShuffleIdents<-function(seu, ident.name, n, random.seed=42, show.progress=TRUE){
    Idents(seu)<-ident.name
    true.ident <- seu@meta.data[ , ident.name, drop=TRUE] # get a vector of cell annotations
    
    set.seed(random.seed)
    if(show.progress==TRUE){
        library(progress)
        pb <- progress_bar$new(total = n,
                               format= "Computing the background [:bar] :percent in :elapsed eta: :eta",
                               clear = TRUE, width= 60,show_after=5)
        for (i in 1:n){
            pb$tick()
            rand.ident.name <- paste("random",i, sep="_")
            seu@meta.data[, rand.ident.name] <- sample(true.ident)
        }
    }
    
    return(seu)
}

LRDB<-read.csv("/Users/chensijie/Codes/DEPair/LR_DASkelly.csv")
#LRDB<-LRDB[1:100,]
load("/Users/chensijie/Codes/DEPair/Colon_0.6_0730.RDS")
immune.combined<-subset(immune.combined, 
                        subset = `celltype0627`!="16.Epithelial cells"&
                          `celltype0627`!="17.Endothelial cells" &
                          `celltype0627`!="18.Fibroblasts")
ls()

LRDB            <- TrimDB.primary(immune.combined,LRDB)
immune.combined <- TrimSeuratFeatures.primary(immune.combined,LRDB)
crosstalk.df    <- GetCrosstalks(seu=immune.combined, LRDB=LRDB, 
                               ident.name="celltype0627",expr.cutoff = 1)
LRDB            <- TrimDB.secondary(crosstalk.df, LRDB)
immune.combined <- TrimSeuratFeatures.secondary(immune.combined, crosstalk.df)


immune.combined <- ShuffleIdents(seu = immune.combined, 
                                ident.name = "celltype0627",n=1000, show.progress = T)
DimPlot(immune.combined, group.by = "celltype0627")


aa=GetRandExpression(seu=immune.combined,gene= "Apoe",rand.batch=1000, celltype="0.IgD+ mature B")
bb=GetRandExpression(seu=immune.combined,gene= "Sorl1",rand.batch=1000, celltype="0.IgD+ mature B")
aa*bb


CountExtremes <- function(rand.batch, seu, celltype.ligand, celltype.receptor, ligand, receptor){
    rand.ligand  <-GetRandExpression(rand.batch, seu=seu,gene=ligand,  celltype=celltype.ligand)
    rand.receptor<-GetRandExpression(rand.batch, seu=seu,gene=receptor,celltype=celltype.receptor)
    rand.intensity<-rand.ligand*rand.receptor
    return(rand.intensity)
}

GetNullDistribution<-function(seu, nPerm, nCores, 
                              celltype.ligand, celltype.receptor, ligand,receptor){
    if (nCores==1){
        background<- lapply(1:nPerm, CountExtremes, seu=seu, 
                          celltype.ligand=celltype.ligand, celltype.receptor=celltype.receptor, 
                          ligand=ligand,receptor=receptor)
        background<- unlist(background)
        return(background)
    }else{
        library(parallel)
        cl <- makeCluster(getOption("cl.cores", nCores), type="FORK")
        clusterExport(cl=cl, varlist=c("seu",
                                       "celltype.ligand",
                                       "celltype.receptor",
                                       "ligand","receptor"), envir=environment())
        
        background<-clusterApply(cl,1:nPerm,
                                 CountExtremes, seu=seu, 
                                 celltype.ligand=celltype.ligand, 
                                 celltype.receptor=celltype.receptor, 
                                 ligand=ligand,receptor=receptor) 
        stopCluster(cl)
        background<- unlist(background)
        return(background)
    }
  
}


ComputePValue<- function(seu, nPerm, nCores, value, 
                         celltype.ligand, celltype.receptor, ligand,receptor,
                         return.background=FALSE){
    bkg<-GetNullDistribution(seu = seu, nPerm = nPerm, nCores = nCores,
                             celltype.ligand =celltype.ligand, 
                             celltype.receptor=celltype.receptor,
                             ligand =ligand,receptor=receptor)
    r <- sum(bkg>=value)
    # North BV, Curtis D, Sham PC. A note on the calculation of empirical P
    # values from Monte Carlo procedures. Am J Hum Genet. 2002;71(2):439-441. 
    # doi:10.1086/341527
    p <- (r+1.0)/(nPerm+1.0)
    if (return.background==FALSE){
        return(p)
    }else{
        return(list(p=p, bkg=bkg))
    }
}

tic("parallel")
GetNullDistribution(seu = immune.combined, nPerm = 400, nCores = 4,
                    celltype.ligand ="0.IgD+ mature B", celltype.receptor="0.IgD+ mature B",
                    ligand = "Apoe",receptor = "Sorl1")
toc()

tic("sequential")
GetNullDistribution(seu = immune.combined, nPerm = 400, nCores = 1,
                    celltype.ligand ="0.IgD+ mature B", celltype.receptor="0.IgD+ mature B",
                    ligand = "Apoe",receptor = "Sorl1")
toc()

ComputePValue(seu = immune.combined, nPerm = 20, nCores = 4, value = 7.375349,
              celltype.ligand ="0.IgD+ mature B",
              celltype.receptor="0.IgD+ mature B",
              ligand = "Apoe",receptor = "Sorl1",
              return.background = TRUE
              )


celltype.interactions<-crosstalk.df %>% group_by(celltype.1,celltype.2) %>% summarise(interaction=sum(intensity))
hist(celltype.interactions$interaction,breaks = 100)
