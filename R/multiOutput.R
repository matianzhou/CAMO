##' Analysis results for multiple pairs: visualization outputs including overall 
##' pathway clustering and output for each pathway
##' The \code{multiOutput} is function to generate visualization outputs 
##' for multiple pairs: including overall pathway clustering outputs, model MDS plot, 
##' model clustering output, heatmap of gene posterior mean, kegg pathway topology for 
##' each pathway
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

multiOutput <- function(mcmc.merge.list,dataset.names,select.pathway.list,ARS_pathway,
                output=c("clustPathway","mdsModel","clustModel","genePM","keggView"),
                hashtb=NULL,pathways=NULL,keggViewSelect = c(1,2),optK=NULL,
                kegg_pathname=NULL,hs_gene_id=NULL) {

### Multiple pairs output including the following outputs:
## pathway clustering (including consensus clust, heatmap, mds, pathway text mining)
## & pathway level: per pathway mdsModel plot, clustModel heatmap, gene heatmap, 
## pathview (KEGG only)
  
  if(length(output)==0 || is.null(output)){
      stop("at least one type of output has to be chosen")  
  }
  
  orig.path <- getwd()
  pathway.name <- names(select.pathway.list)
  K <- length(pathway.name)
  
  ARS.mat <- ARS_pathway[["ARS.mat"]]
  ARSpvalue.mat <- ARS_pathway[["ARSpvalue.mat"]]
  
  M <- length(mcmc.merge.list)
  P <- ncol(ARS.mat)
  
  if("clustPathway" %in% output){
    dir.path <- "clustPathway"
    if (!file.exists(dir.path)) dir.create(dir.path)
    setwd(paste(orig.path,"/",dir.path,sep=""))
    
    #1. consensus clustering
    results <- clustPathway(ARSpvalue.mat)
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
    res <- mdsPathway(arsPvalue=ARSpvalue.mat,
                          cluster.assign=cluster.assign,
                          scatter.index=scatter.index)
    #5. text mining
    if(!is.null(hashtb) && !is.null(pathways)){
      textMine(hashtb=hashtb,pathways=pathways,cluster.assign=cluster.assign)
    } ## should exclude the scattered pathways
    
    #6. heatmap
    
    res <- heatmapPathway(arsPvalue=ARSpvalue.mat,
                          cluster.assign=cluster.assign,
                          scatter.index=scatter.index)
  }
  setwd(orig.path)
  
  if("mdsModel" %in% output){
    dir.path <- "mdsModel"
    if (!file.exists(dir.path)) dir.create(dir.path)
    setwd(paste(orig.path,"/",dir.path,sep=""))
    
    for(k in 1:K){
      print(paste("mdsModel",k,sep=":"))
      pathk.name <- pathway.name[k]
      res <- mdsModel(unlist(c(ARS.mat[k,])),dataset.names,pathk.name,sep="_")
    }  
  }    
  setwd(orig.path)
  
  if("clustModel" %in% output){
    dir.path <- "clustModel"
    if (!file.exists(dir.path)) dir.create(dir.path)
    setwd(paste(orig.path,"/",dir.path,sep=""))
    
    for(k in 1:K){
      print(paste("clustModel",k,sep=":"))
      pathk.name <- pathway.name[k]
      cluster.assign <- try(SA_algo(unlist(c(ARSpvalue.mat[k,])),dataset.names,sep="_"))
      if(length(unique(cluster.assign))>1 && class(cluster.assign) != "try-error" ){
      res <- clustModel(unlist(c(ARSpvalue.mat[k,])),dataset.names,cluster.assign,
                        pathk.name,sep="_")
      }
      if(length(unique(cluster.assign))==1 || class(cluster.assign) == "try-error" ){
        warning("clustModel only identifies one cluster")
      res <- clustModelOne(unlist(c(ARSpvalue.mat[k,])),dataset.names,
                          pathk.name,sep="_") 
      }
    }  
  }
  setwd(orig.path)
  
  if("genePM" %in% output){
    dir.path <- "genePM"
    if (!file.exists(dir.path)) dir.create(dir.path)
    setwd(paste(orig.path,"/",dir.path,sep=""))
    for(k in 1:K){
      print(paste("genePM",k,sep=":"))
      pathk.name <- pathway.name[k]
      pathway.genes <- select.pathway.list[[k]]
      signPM.list <- lapply(mcmc.merge.list,function(x) apply(x,1,mean))
      names(signPM.list) <- dataset.names
      hm <- genePM(signPM.list, pathway.genes=pathway.genes, 
                        pathway.name=pathk.name)
    }
  }
  setwd(orig.path)
  
  if("keggView" %in% output){
    if(sum(grepl("KEGG",pathway.name))==0) {
      warning("No KEGG pathways") 
    } else{
    dir.path <- "keggView"
    if (!file.exists(dir.path)) dir.create(dir.path)
    setwd(paste(orig.path,"/",dir.path,sep=""))
    
    kegg.pathway.name <- pathway.name[grep("KEGG",pathway.name)]
    K_KEGG <- length(kegg.pathway.name)
    dat1 <- mcmc.merge.list[[keggViewSelect[1]]]
    dat2 <- mcmc.merge.list[[keggViewSelect[2]]]
    for(k in 1:K_KEGG){
      print(paste("keggView",k,sep=":"))
      keggk.name <- kegg.pathway.name[k]
      overlap.genes <- intersect(rownames(dat1),select.pathway.list[[keggk.name]])
      signPM.mat <- cbind(apply(dat1[overlap.genes,],1,mean),
                          apply(dat2[overlap.genes,],1,mean))
      colnames(signPM.mat) <- dataset.names[keggViewSelect]
      keggk.name1 <- gsub("KEGG ","",keggk.name) 
      #library(KEGG.db)
      #xx <- unlist(as.list(KEGGPATHID2NAME))
      pathwayID <- names(kegg_pathname)[which(kegg_pathname==keggk.name1)]
      res <- keggView(mat=signPM.mat,pathwayID)
      if(grepl("/",keggk.name)){
        keggk.name <- gsub("/","-",keggk.name)
      }
      hsaName <- paste("hsa",pathwayID,sep="")
      file.rename(paste(hsaName,"..multi.png",sep=""), 
                  paste(keggk.name,".png",sep=""))
      file.remove(paste(hsaName,".xml",sep=""))
      file.remove(paste(hsaName,".png",sep=""))
     }
   } 
  }  
    setwd(orig.path)       
    
  print("Multiple pairs analysis completed.") 
  
}


clustPathway <- function(arsPvalue) {
  #require(ConsensusClusterPlus)
  #set your working dir, automatically save there
  
  results = ConsensusClusterPlus(d=t(-log10(arsPvalue)),maxK=10,reps=50,pItem=0.8,pFeature=1,title="Pathway clustering",clusterAlg="hc",innerLinkage="ward.D2",finalLinkage="ward.D2",seed=15213,plot="png")
  
  return(results) ## a list of K elements (each k represents number of clusters)
}

clustNumber <- function(results){
  Kvec = 2:length(results);
  x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
  PAC = rep(NA,length(Kvec)) 
  names(PAC) = paste("K=",Kvec,sep="") # from 2 to 10 (maxK)
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i
# The optimal K
  optK = Kvec[which.min(PAC)-2]
  return(optK)
}

mdsPathway <- function(arsPvalue,cluster.assign,scatter.index=NULL) {
  
  ## plot MDS for all pathways
  ## arsPvalue is a matrix of K (pathways) rows and choose(M,2) columns 
  ## cluster.assign is the result from consensus clustering
  
  C <- length(unique(cluster.assign)) #number of clusters
  if(C>9){
    warning("Too many clusters, not enough colors")
  }
  
  dist.mat <-  dist(-log10(arsPvalue),method = "euclidean",
                    upper = TRUE, diag = TRUE)
  
  fit <- cmdscale(dist.mat,k=2)
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
  #pdf(paste("mdsPathway","_K_",C,".pdf",sep=""))
  jpeg(paste("mdsPathway","_K_",C,".jpeg",sep=""),quality = 100)
  
  p <- ggplot() +
       ggtitle("") +
    xlab("Coordinate 1") + ylab("Coordinate 2") + 
    xlim(c(-xlimit,xlimit)) + ylim(c(-ylimit,ylimit)) + 
    geom_point(aes(x, y), shape=shapes, 
               color = colors ,size=sizes) +
    geom_point(aes(xcenter,ycenter),
               shape=unique.shape, color = unique.color, 
               size =5) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 15, hjust=0.5,face="bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  print(p)
  dev.off()
  return(p)
}  

textMine <- function(hashtb,pathways,cluster.assign){
  ##cluster.assign with pathway names (w/o scatterness)
  
  tmk <- TextMine(hashtb=hashtb, pathways= pathways,
                  pathway=names(cluster.assign), result=cluster.assign)
  C <- length(unique(cluster.assign))
  tm_filtered <- list()
  for (i in 1:C){ 
    tm_filtered[[i]] <- tmk[[i]][which((as.numeric(tmk[[i]][,4]) < 0.05)), ]
  }
  
  pathway.summary <- lapply(1:C, function(x) names(which(cluster.assign==x)))
  writeTextOut(tm_filtered,C,pathway.summary)
}


heatmapPathway <- function(arsPvalue, cluster.assign,scatter.index=NULL){
  
  ## cluster.assign is the result from consensus clustering
  arsPvalue <- data.matrix(arsPvalue)
  C <- length(unique(cluster.assign)) #number of clusters
  dataOrder <- -log10(arsPvalue)[unlist(sapply(1:C,function(x) which(cluster.assign==x))),]
  colnames(dataOrder) <- colnames(arsPvalue)
  row.sep <-  c(0,cumsum(unlist(sapply(1:C,function(x) sum(cluster.assign==x)))) )   
  
  if(!is.null(scatter.index)){
    cluster.assign[scatter.index] <- -1
    dataOrder <- -log10(arsPvalue)[unlist(sapply(c(1:C,-1),function(x) which(cluster.assign==x))),]      
    colnames(dataOrder) <- colnames(arsPvalue)
    row.sep <-  c(0,cumsum(unlist(sapply(c(1:C,-1),function(x) sum(cluster.assign==x)))) )                    
  }
  #rownames(dataOrder) <- sapply(rownames(dataOrder),function(x) substr(x,1,15))
  ordered.cluster.assign <- cluster.assign[rownames(dataOrder)]
  row.colors <- rep(NA, length(ordered.cluster.assign) )
  for(i in 1:length(row.colors)){
    if(ordered.cluster.assign[i]== -1){
      row.colors[i] <- "gray"
    } else {
      row.colors[i] <- rainbow(C)[ordered.cluster.assign[i]]
    }
  }
  #pdf(paste("heatmapPathway","_K_",C,".pdf",sep=""))
  jpeg(paste("heatmapPathway","_K_",C,".jpeg",sep=""),quality = 100)
  par(cex.main=1, font.lab=2, font.axis=2)
  hm<-heatmap.2(dataOrder, symm=F,main=NULL,
                cexCol=0.7,cexRow=0.3,adjCol= c(NA,-1),
                rowsep=row.sep,
                sepwidth=c(0.1, 0.3),  # width of the borders
                sepcolor=c('white'),scale='none',
                symbreaks=T,key=T, keysize=1,symkey=F, 
                dendrogram=c('none'),density.info="none", 
                trace="none",Rowv=F,Colv=T,
                srtCol=50,RowSideColors=row.colors,
                col=greenred,breaks=seq(0,max(dataOrder),by=0.01),
                key.ytickfun=function(){
                  side = 2
                } )
  dev.off()
  
  return(hm)
  
}  

