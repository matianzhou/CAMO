##' Plot mds plot (function only for Rshiny)
##' @title Analysis results for multiple pairs: visualization outputs  
##' @param mcmc.merge.list: a list of merged MCMC output matrices.
##' @param dataset.names: a vector of dataset names.
##' @param select.pathway.list: a list of selected pathways (containing gene 
##' components).
##' @param ARS_pathway: a list of two data frames: pathway specific ARS values and 
##' their permuted p-value (pathway on rows, column being ARS value or the p-values).
##' @param output: five options: "clustPathway" (pathway clustering),"mdsModel"(model 
##' MDS plot),"clustModel" (model clustering output), "genePM" (generating heatmap of 
##' gene posterior mean),"keggView" (generating kegg pathway topology, human KEGG 
##' only). For details, please refer to manuscript. cannot be empty.
##' @param hashtb: hash table for text mining.
##' @param pathways: complete pathway names for text mining.
##' @param keggViewSelect: which two datasets to view in KEGG topology. 
##' @param optK: Optimal number of clusters based on clustering diagnostic results. For
##' "clustPathway" output only. 
##' @param kegg_pathname: KEGG pathway name list. For "keggView" only.
##' @param hs_gene_id: Human sapiens gene id. For "keggView" only.

##' @return stored output in created folders.
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step
##' #select.pathway.list from the pathSelect step
##' #ARS_pathway from the multiARS step
##' data(hashtb) #include hashtb & pathways
##' dataset.names <- c("hb","hs","ht","mb","ms","mt")
##' library(KEGG.db)
##' kegg_pathname <- unlist(as.list(KEGGPATHID2NAME))
##' library("org.Hs.eg.db")
##' hs_gene_id <- unlist(mget(x=rownames(mcmc.merge.list[[1]]),
##' envir=org.Hs.egALIAS2EG))
##' multiOutput(mcmc.merge.list,dataset.names,select.pathway.list,
##' ARS_pathway, output=c("clustPathway","mdsModel","clustModel","genePM","keggView"),
##' hashtb=hashtb,pathways=pathways,keggViewSelect = c(1,4),optK=7)
##' }

mdsPlotRShiny <- function(mcmc.merge.list, dataset.names, select.pathway.list,
  ARS_pathway, hashtb=NULL, pathways=NULL,  optK=NULL) {
  
  pathway.name <- names(select.pathway.list)
  K <- length(pathway.name)
  
  ARS.mat <- ARS_pathway[["ARS.mat"]]
  ARSpvalue.mat <- ARS_pathway[["ARSpvalue.mat"]]
  
  M <- length(mcmc.merge.list)
  P <- ncol(ARS.mat)
      
  #1. consensus clustering
  results <- ConsensusClusterPlus(d=t(-log10(ARSpvalue.mat)),maxK=10,reps=50,pItem=0.8,pFeature=1,title="Pathway clustering",clusterAlg="hc",innerLinkage="ward.D2",finalLinkage="ward.D2",seed=15213,plot="png")
    
  #2. determine optimal number of clusters
  if(is.null(optK)){
    optK <- clustNumber(results)
  }
  cluster.assign <- results[[optK]]$consensusClass

  #3. identify scattered objects
  scatter.index <- scatter(-log10(ARSpvalue.mat), cluster.assign)
  if(length(scatter.index) >= round(K/5)){
    scatter.index <- NULL
  }

  #4. mds plot
  C <- length(unique(cluster.assign)) #number of clusters
  if(C>9){
    warning("Too many clusters, not enough colors")
  }
  
  dist.mat <-  dist(-log10(ARSpvalue.mat),method = "euclidean",
                    upper = TRUE, diag = TRUE)
  fit <- cmdscale(dist.mat,k=2)
  coordi <- fit

  x <- fit[,1]
  y <- fit[,2]
  xlimit <- ifelse(abs(min(x))>abs(max(x)),abs(min(x)),abs(max(x)))
  ylimit <- ifelse(abs(min(y))>abs(max(y)),abs(min(y)),abs(max(y)))
  xcenter <- tapply(x,as.factor(cluster.assign),mean)
  ycenter <- tapply(y,as.factor(cluster.assign),mean)
  
  if(!is.null(scatter.index)){
    cluster.assign[scatter.index] <- -1
  }
  
  unique.color <- rainbow(C)
  unique.shape <- 1:C
  sizes <- shapes <- colors <- cluster.assign
  for(i in 1:(C+1)){
    colors[cluster.assign==i] <- unique.color[i]
    shapes[cluster.assign==i] <- unique.shape[i]
    sizes[cluster.assign==i] <- 2
    if(i== (C+1)){
      colors[cluster.assign== -1] <- "gray50"
      shapes[cluster.assign== -1] <- 20 
      sizes[cluster.assign== -1] <- 1
    }
  }
  
  #5. text mining
  if(!is.null(hashtb) && !is.null(pathways)){
    tmk <- TextMine(hashtb=hashtb, pathways= pathways,
                    pathway=names(cluster.assign), result=cluster.assign)
    C <- length(unique(cluster.assign))
    tm_filtered <- list()
    for (i in 1:C){ 
      tm_filtered[[i]] <- tmk[[i]][which((as.numeric(tmk[[i]][,4]) < 0.05)), ]
    }
  }
  
   return(list(Xlim=xlimit, Ylim=ylimit, Shapes=shapes, Colors=colors, 
    Sizes=sizes, Xcenter=xcenter, Ycenter=ycenter, Unique.shape=unique.shape, 
    Unique.color=unique.color, 
    Cluster.assign=cluster.assign, Coordi=coordi, 
    Tm_filtered=tm_filtered)) 
}


