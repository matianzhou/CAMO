##' Merge multiple MCMCout datasets
##' The \code{mergeMCMC} is function to merge multiple MCMCout datasets by
##' matching orthologs
##' @title Merge multiple MCMCout datasets preparing for cross-species
##' analysis
##' @param mcmc.list: a list of MCMC output matrices.
##' @param species: a vector specie names of same length as mcmc.list.
##' @param ortholog.db: the ortholog object (in R environment)
##' @param ortholog.file: the ortholog file to be imported
##' @param reference: the index of the reference data, the outputted merged 
##' list will be named using the rownames of this data.

##' @return an merged list of multiple MCMCout datasets (with same number of 
##' rows and rownames)
##' @export
##' @examples
##' \dontrun{
##' data(hb)
##' summaryDE <- indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' hb_pData <- summaryDE[,c(3,1)]
##' hb_MCMCout <- bayes(hb_pData, seed=12345)
##' data(hs)
##' summaryDE <- indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' hs_pData <- summaryDE[,c(3,1)]
##' hs_MCMCout <- bayes(hs_pData, seed=12345)
##' data(ht)
##' summaryDE <- indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' ht_pData <- summaryDE[,c(3,1)]
##' ht_MCMCout <- bayes(ht_pData, seed=12345)
##' data(mb)
##' summaryDE <- indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' mb_pData <- summaryDE[,c(3,1)]
##' mb_MCMCout <- bayes(mb_pData, seed=12345)
##' data(ms)
##' summaryDE <- indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' ms_pData <- summaryDE[,c(3,1)]
##' ms_MCMCout <- bayes(ms_pData, seed=12345)
##' data(mt)
##' summaryDE <- indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' mt_pData <- summaryDE[,c(3,1)]
##' mt_MCMCout <- bayes(mt_pData, seed=12345)
##' 
##' #1. single pair example
##' mcmc.list <- list(hb_MCMCout,mb_MCMCout)
##' species <- c("human","mouse")
##' data(hm_orth)
##' mcmc.merge.list <- mergeMCMC(mcmc.list,species = species,
##' ortholog.db = hm_orth, reference=1) 
##' 
##' #2. multiple pairs example
##' mcmc.list <- list(hb_MCMCout,hs_MCMCout,ht_MCMCout,
##'                   mb_MCMCout,ms_MCMCout,mt_MCMCout)
##' species <- c(rep("human",3), rep("mouse",3))
##' data(hm_orth)
##' mcmc.merge.list <- mergeMCMC(mcmc.list,species = species,
##' ortholog.db = hm_orth, reference=1)                
##' }

mergeMCMC <- function(mcmc.list,species,ortholog.db,ortholog.file=NULL,
                  reference=1){

  M <- length(mcmc.list)
  mcmc.merge.list <- DEgene.merge.list <- vector("list",M)
  
  if(is.null(ortholog.file)){
    ortholog.data <- ortholog.db
  } else {
    ortholog.data <- read.csv(ortholog.file,stringsAsFactors = F,header=T) 
  }
    
       gene.list <- lapply(mcmc.list,rownames)
       DEgene.list <- lapply(mcmc.list,function(x) rownames(x)[attr(x,"DEindex")])

       match_gene <- orthMatch(gene.list,species,ortholog.data)
       ref_gene <- match_gene[[reference]] 
       ## match_gene and ref_gene of same dimension
       for(m in 1:M){
         mcmc.merge.list[[m]] <- mcmc.list[[m]][match_gene[[m]],]
         rownames(mcmc.merge.list[[m]]) <- ref_gene
         DEgenes <- ref_gene[which(match_gene[[m]] %in% DEgene.list[[m]])]
         DEindex <- which(rownames(mcmc.merge.list[[m]])%in% DEgenes)
         attr(mcmc.merge.list[[m]],"DEindex") <- DEindex
       }

   return(mcmc.merge.list)
}


orthMatch <- function(gene.list,species,ortholog.data){
  M <- length(gene.list)
  index.out <- gene.out <- vector("list",M)
  for(m in 1:M){
    spec <- species[m]
    gene <- gene.list[[m]]
    index.out[[m]] <- which(ortholog.data[,spec] %in% gene)
  }
    common.index <- Reduce(intersect,index.out)
  
  for(m in 1:M){
    spec <- species[m]
    gene.out[[m]] <- ortholog.data[common.index,spec]
  }

  return(gene.out)
}