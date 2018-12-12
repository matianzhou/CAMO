##' Merge multiple raw datasets
##' The \code{mergeRaw} is function to merge multiple raw datasets by
##' matching orthologs
##' @title Merge multiple raw datasets preparing for cross-species
##' analysis
##' @param data.list: a list of data matrices.
##' @param species: a vector specie names of same length as data.list.
##' @param ortholog.db: the ortholog object (in R environment)
##' @param ortholog.file: the ortholog file to be imported
##' @param reference: the index of the reference data, the outputted merged 
##' list will be named using the rownames of this data.
##' @param unique: logical value indicating whether to take only the "one to one" unique match (T) or allow "multiple to one" and "one to multiple" match (F) 
##' @return an merged list of rdatasets (with same number of 
##' rows and rownames)
##' @export
##' @examples
##' \dontrun{
##' data(hb)
##' hb.data <- data
##' hb.group <- group
##' data(hs)
##' hs.data <- data
##' hs.group <- group
##' data(ht)
##' ht.data <- data
##' ht.group <- group
##' data(mb)
##' mb.data <- data
##' mb.group <- group
##' data(ms)
##' ms.data <- data
##' ms.group <- group
##' data(mt)
##' mt.data <- data
##' mt.group <- group
##' data.list <- list(hb.data,hs.data,ht.data,mb.data,ms.data,mt.data)
##' species <- c(rep("human",3), rep("mouse",3))
##' data(hm_orth)
##' data.merge.list <- mergeRaw(data.list,species=species,
##' ortholog.db = hm_orth, reference=1,unique=T)
##' }


mergeRaw <- function(data.list,species,ortholog.db,ortholog.file=NULL,
                     reference=1,unique=T){
  
  M <- length(data.list)
  data.merge.list <- vector("list",M)
  gene.list <- lapply(data.list,rownames)
  
  if(is.null(ortholog.file)){
    ortholog.data <- ortholog.db
  } else {
    ortholog.data <- read.csv(ortholog.file,stringsAsFactors = F,header=T) 
  }
  unique.index <- which(!duplicated(apply(ortholog.data,1,paste,collapse=",")))
  ortholog.data <- ortholog.data[unique.index,]
  
  if(unique){
    ref.spec <- species[[reference]]
    ref.orth <- ortholog.data[,ref.spec]
    ref.gene <- intersect(gene.list[[reference]],names(which(table(ref.orth)==1)))
    spec.gene <- vector("list",M)
    
    ## take common first, then match orthologs
    
    for(m in 1:M){
      spec <- species[m]
      spec.gene[[m]] <- ortholog.data[which(ortholog.data[,spec] %in% rownames(data.list[[m]])),ref.spec]
    }
    common.gene <- intersect(ref.gene,Reduce(intersect,spec.gene))
    
    data.new.list <- lapply(1:M,function(m) {
      a <- matrix(NA,nrow=length(common.gene),
                  ncol=ncol(data.list[[m]]));
      rownames(a) <- common.gene;
      colnames(a) <- colnames(data.list[[m]]);
      return(a) } )
    
    m_ref <- which(species==ref.spec)
    m_notref <- which(species!=ref.spec)
    
    for(m in m_ref){
      data.new.list[[m]] <- data.list[[m]][common.gene,]
    }
      
    for(m in m_notref){
      spec <- species[m]
      index <- which(ortholog.data[,ref.spec] %in% common.gene)
      names(index) <- ortholog.data[index,ref.spec]
      order.index <- index[common.gene]
      data.new.list[[m]][common.gene,] <- data.list[[m]][ortholog.data[order.index,spec],]
    }  
    
    data.merge.list <- data.new.list
  }
  
  if(!unique){
    
    ref.spec <- species[[reference]]
    ref.orth <- ortholog.data[,ref.spec]
    ref.gene <- intersect(gene.list[[reference]],ref.orth)
    unique.ref.gene <- intersect(gene.list[[reference]],names(which(table(ref.orth)==1)))
    spec.gene <- vector("list",M)
    
    ## take common first, then match orthologs
    
    for(m in 1:M){
      spec <- species[m]
      spec.gene[[m]] <- ortholog.data[which(ortholog.data[,spec] %in% rownames(data.list[[m]])),ref.spec]
    }
    common.gene <- intersect(ref.gene,Reduce(intersect,spec.gene))
    
    data.new.list <- lapply(1:M,function(m) {
      a <- matrix(NA,nrow=length(common.gene),
                  ncol=ncol(data.list[[m]]));
      rownames(a) <- common.gene;
      colnames(a) <- colnames(data.list[[m]]);
      return(a) } )
    
    m_ref <- which(species==ref.spec)
    m_notref <- which(species!=ref.spec)
    
    for(m in m_ref){
      data.new.list[[m]] <- data.list[[m]][common.gene,]
    }
    
    g1.index <- which(common.gene %in% unique.ref.gene)
    g2.index <- (1:length(common.gene))[-g1.index]
    
    for(m in m_notref){
      spec <- species[m]
      index <- which(ortholog.data[,ref.spec] %in% common.gene[g1.index])
      names(index) <- ortholog.data[index,ref.spec]
      order.index <- index[common.gene[g1.index]]
      data.new.list[[m]][common.gene[g1.index],] <- data.list[[m]][ortholog.data[order.index,spec],]
    }  
    
    for(g in g2.index){
      gene <- common.gene[g]
      row <- which(ref.orth==gene)
        for(m in m_notref){
          spec <- species[m]
          spec.gene <- ortholog.data[row,spec]
          match.gene <- intersect(rownames(data.list[[m]]),spec.gene)
          if(length(match.gene)>1){
            iqr <- apply(data.list[[m]][match.gene,],1,IQR)
            data.new.list[[m]][g,] <- data.list[[m]][match.gene[which.max(iqr)],]
          } else if(length(match.gene)==1){
            data.new.list[[m]][g,] <- data.list[[m]][match.gene,]
          }
        }
    }

      data.merge.list <- data.new.list
    
    }
    
  return(data.merge.list)
}

