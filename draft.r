ShuffleIdent<- function(random.seed, seurat.obj, ident.name, expr.cutoff, lig.name, rec.name){
    # set the random seed
    set.seed(random.seed) 
    
    
    # shuffle the identities
    Idents(seurat.obj)<-ident.name
    identity <- seurat.obj@meta.data[ , ident.name, drop=TRUE] # get a vector of cell annotations
    seurat.obj$random<-sample(identity)      # randomly permute the annotations
    
    
    # compute the average expression after permutation
    Idents(seurat.obj) <- "random"
    expr <- AverageExpression(object = seu, 
                              assays = "RNA" ,slot = "data", verbose = F)$RNA # avg expr
    
    
    
    expr[expr<expr.cutoff] <- 0              # drop low expressions
    expr <- as(as.matrix(expr), "dgCMatrix") # sparsify the matrix
    
    
    result<-list("expr.ligand"  =expr[lig.name,  ], 
                 "expr.receptor"=expr[rec.name,  ]) 
    return(result)
}

ComputeBackground <- function(seu, n, LRDB, group.by, expr.cutoff=1, cores=4){
    
    message("Triming LR-pair DB")
    # trim ligand-receptor pair database 
    ligand.indata   <- unique(LRDB$ligand_symbol[LRDB$ligand_symbol %in%rownames(seu)])
    receptor.indata <- unique(LRDB$receptor_symbol[LRDB$receptor_symbol %in%rownames(seu)])
    lrpair.indata   <-(LRDB$ligand_symbol %in% ligand.indata)&(LRDB$receptor_symbol %in% receptor.indata)
    LRDB            <- LRDB[lrpair.indata, ]

    message("Triming input object")
    # trim seurat object
    features.sel<-unique(c(LRDB$ligand_symbol, LRDB$receptor_symbol))
    seu<-subset(seu, features = unique(c(LRDB$ligand_symbol, LRDB$receptor_symbol)))
    for (assay_name in names(seu@assays)){
        if (assay_name!="RNA"){ seu@assays[[assay_name]]<-NULL }
    }
    Idents(seu)<-group.by
    
    if (cores>1){
        message("Computing parallely")
        library(parallel)
        cl <- makeCluster(getOption("cl.cores", cores), type="FORK")
        clusterExport(cl=cl, varlist=c("seu","n","LRDB","group.by", "expr.cutoff"), envir=environment())
        Backgrounds<-clusterApply(cl,1:n,ShuffleIdent, 
                                    seurat.obj  = seu, 
                                    ident.name  = group.by,
                                    expr.cutoff = expr.cutoff, 
                                    lig.name    = LRDB$ligand_symbol, 
                                    rec.name    = LRDB$receptor_symbol)   #加减乘除
        stopCluster(cl)
        return (Backgrounds)
    }else{
        message("Computing Sequentially")
        NPerm <- n
        Backgrounds<-list()
        library(progress)
        pb <- progress_bar$new(total = NPerm,
                               format= "Computing the background [:bar] :percent in :elapsed eta: :eta",
                               clear = T, width= 60,show_after=50)
        for (i in 1:n){
            pb$tick()
            Backgrounds[[i]]<-ShuffleIdent(random.seed = i, 
                                           seurat.obj  = seu, 
                                           ident.name  = group.by,
                                           expr.cutoff = expr.cutoff, 
                                           lig.name    = LRDB$ligand_symbol,
                                           rec.name    = LRDB$receptor_symbol)
        }
        return (Backgrounds)
    }
    

}







PairMatch <-function(expr.ligand, expr.receptor, celltype.1, celltype.2, LRDB){
    message(paste(celltype.1, "-> ⊃-", celltype.2))
    inner.prod <- expr.ligand[,celltype.1,drop=FALSE] * expr.receptor[,celltype.2,drop=FALSE]
    row.names(inner.prod) <- LRDB$pair
    pos.sel    <- inner.prod[,1]>0
    if (sum(pos.sel)==0){ 
        return (NULL)
    }
    
    if (sum(pos.sel)>1){ 
        positive.scores <- inner.prod[pos.sel,,drop=TRUE]
        LR.scores  <- sort(positive.scores, decreasing = T)   
    }else{ #一个数是一个退化的vector，默认没有name，单行df在drop后不形成vector，变成一个无name的数，故此特判
        positive.scores <- inner.prod[pos.sel,,drop=FALSE]
        LR.scores       <- positive.scores[1,1]
        names(LR.scores)<-rownames(positive.scores) 
    }
    return(LR.scores)
}


CountExtremeEvent<-function(shuffle, LR.scores, 
                            celltype.1, celltype.2, 
                            expr.ligand, expr.receptor, LRDB){
    # compute the interaction scores in shuffled data
    expr.ligand    <- shuffle$expr.ligand
    expr.receptor  <- shuffle$expr.receptor
    
    inner.prod     <- expr.ligand[,celltype.1,drop=FALSE] * expr.receptor[,celltype.2,drop=FALSE]
    row.names(inner.prod) <- LRDB$pair
    
    
    LR.scores.null <- inner.prod[names(LR.scores),]
    extreme        <- as.numeric(LR.scores.null>LR.scores)
    return(extreme)
}


ComputePVal <- function(LR.scores, Backgrounds, 
                        celltype.1, celltype.2, 
                        expr.ligand, expr.receptor, LRDB, cores=1){
    n <-length(Backgrounds)
    
    if (cores>1){
        library(parallel)
        cl <- makeCluster(getOption("cl.cores", cores), type="FORK")
        clusterExport(cl=cl, varlist=c("LR.scores","Backgrounds","celltype.1","celltype.2",
                                       "expr.ligand","expr.receptor","LRDB"), envir=environment())
        
        count.list<-clusterApply(cl,Backgrounds,
                           CountExtremeEvent, 
                             LR.scores=LR.scores, celltype.1=celltype.1, celltype.2=celltype.2, 
                             expr.ligand=expr.ligand, expr.receptor=expr.receptor, LRDB=LRDB) 
        stopCluster(cl)
    }else{
        count.list <-lapply(Backgrounds, CountExtremeEvent, 
                                     LR.scores=LR.scores, celltype.1=celltype.1, celltype.2=celltype.2, 
                                     expr.ligand=expr.ligand, expr.receptor=expr.receptor, LRDB=LRDB)    
    }

    LR.scores.pval <- Reduce("+",count.list) / n
    return(LR.scores.pval)
}


ComputeCrosstalkScores <- function(seu, LRDB, group.by, Backgrounds, expr.cutoff=1, cores=1){
    
    # trim ligand-receptor pair database 
    ligand.indata   <- unique(LRDB$ligand_symbol[LRDB$ligand_symbol %in%rownames(seu)])
    receptor.indata <- unique(LRDB$receptor_symbol[LRDB$receptor_symbol %in%rownames(seu)])
    lrpair.indata   <-(LRDB$ligand_symbol %in% ligand.indata)&(LRDB$receptor_symbol %in% receptor.indata)
    LRDB            <- LRDB[lrpair.indata, ]
    

    # trim seurat object
    features.sel<-unique(c(LRDB$ligand_symbol, LRDB$receptor_symbol))
    seu<-subset(seu, features = unique(c(LRDB$ligand_symbol, LRDB$receptor_symbol)))
    for (assay_name in names(seu@assays)){
        if (assay_name!="RNA"){ seu@assays[[assay_name]]<-NULL }
    }
    
    # obtain expression vectors for ligands and corresponding receptors
    Idents(seu)<-group.by
    expr <- AverageExpression(object = seu, assays = "RNA" ,slot = "data",verbose = F)$RNA
    expr[expr<expr.cutoff] <- 0
    expr <- as(as.matrix(expr), "dgCMatrix") 
    expr.ligand  <-expr[LRDB$ligand_symbol,  ]
    expr.receptor<-expr[LRDB$receptor_symbol,]

    library(dplyr)
    crosstalk<-list()
    for (ct1 in levels(seu)){
        for (ct2 in levels(seu)){
            LR.scores <- PairMatch(expr.ligand=expr.ligand, 
                                   expr.receptor=expr.receptor,
                                   celltype.1=ct1, celltype.2=ct2, LRDB=LRDB)
            if (is.null(LR.scores)==FALSE){
                LR.scores.pval <- ComputePVal(LR.scores = LR.scores, Backgrounds = Backgrounds, 
                                          celltype.1=ct1,  celltype.2=ct2,
                                          expr.ligand=expr.ligand, expr.receptor=expr.receptor,
                                          LRDB=LRDB, cores=cores)
                crosstalk[[ct1]][[ct2]] <- data.frame(score=LR.scores, pval=LR.scores.pval)%>%arrange(pval)
            }else{
                crosstalk[[ct1]][[ct2]] <- data.frame(score=numeric(), pval=numeric())
            }
            
            
            
        }
    }
    return(crosstalk)
}
