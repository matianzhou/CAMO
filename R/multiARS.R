##' Resemblance analysis for multiple pairs: global ARS and its permuted p-value
##' The \code{multiARS_global} is function to perform resemblance analysis 
##' for multiple pairs, generating global ARS and its permuted p-value. 
##' @title Resemblance analysis for multiple pairs: global ARS and its 
##' permuted p-value.
##' @param mcmc.merge.list: a list of merged MCMC output matrices.
##' @param dataset.names: a vector of dataset names.
##' @param measure: three types of ARS measures to be used: "youden",
##' "Fmeasure","geo.mean". Default is "Fmeasure". 
##' @param parallel: whether to perform parallel computing in permutation.
##' @param cpu: if parallel=T, how many cpus to be used.
##' @param B: number of permutations. 

##' @return two lists: global ARS values and its permuted p-value, in addition,
##' the two data matrices are written to the folder named "arsGlobal".
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step (example 1)
##' dataset.names <- c("hb","hs","ht","mb","ms","mt")
##' ARS_global <- multiARS_global(mcmc.merge.list,dataset.names,B=100)
##' }

multiARS_global <- function(mcmc.merge.list,dataset.names,
                             measure="Fmeasure",parallel=F,cpu=2,B=50){
  
  names(mcmc.merge.list) <- dataset.names
  M <- length(mcmc.merge.list)
  P <- choose(M,2)
  ARS <- ARSpvalue <- matrix(NA,M,M)
  rownames(ARS) <- colnames(ARS) <- 
  rownames(ARSpvalue) <- colnames(ARSpvalue) <- dataset.names
  diag(ARS) <- diag(ARSpvalue) <- 1
  
  for(i in 1:(M-1)){
    for(j in (i+1):M){
      dat1 <- mcmc.merge.list[[i]]
      dat2 <- mcmc.merge.list[[j]]
      deIndex1 <- attr(dat1,"DEindex") 
      deIndex2 <- attr(dat2,"DEindex") 
      if(parallel==F){
        ARS[i,j] <- ARS[j,i] <- ars(dat1,dat2,deIndex1,deIndex2,measure=measure)
        ARSperm <- permGlobal(dat1,dat2,B=B)
        ARSpvalue[i,j] <- ARSpvalue[j,i] <- arsPvalue(ARS[i,j],ARSperm)
      } else if(parallel==T){
        ARS[i,j] <- ARS[j,i] <- ars(dat1,dat2,deIndex1,deIndex2,measure=measure)        
        ARSperm <- permParGlobal(dat1,dat2,B=B,cpu=cpu)
        ARSpvalue[i,j] <- ARSpvalue[j,i] <- arsPvalue(ARS[i,j],ARSperm)
      }
      print(paste("pair: dataset ",i," and dataset ",j,sep=""))
    }
  }
  dir.path <- "arsGlobal"
  if (!file.exists(dir.path)) dir.create(dir.path)
  write.csv(ARS,file=paste(paste(dir.path,"ARS_global_",sep="/"),M,".csv",sep=""))
  write.csv(ARSpvalue,file=paste(paste(dir.path,"ARSpvalue_global_",sep="/"),M,".csv",sep=""))
  out <- list(ARS=ARS,ARSpvalue=ARSpvalue)
  return(out)
}

##' Resemblance analysis for multiple pairs: pathway specific ARS and their 
##' permuted p-value
##' The \code{multiARS_pathway} is function to perform resemblance analysis 
##' for multiple pairs, generating pathway specific ARS and their permuted 
##' p-value. 
##' @title Resemblance analysis for multiple pairs: pathway specific ARS and 
##' their permuted p-value.
##' @param mcmc.merge.list: a list of merged MCMC output matrices.
##' @param dataset.names: a vector of dataset names.
##' @param select.pathway.list: a list of selected pathways (containing gene 
##' components).
##' @param measure: three types of ARS measures to be used: "youden",
##' "Fmeasure","geo.mean". Default is "Fmeasure".
##' @param parallel: whether to perform parallel computing in permutation.
##' @param cpu: if parallel=T, how many cpus to be used.
##' @param B: number of permutations. 

##' @return a list of two data frames: pathway specific ARS values and their 
##' permuted p-value (pathway on rows, column being ARS value or the 
##' p-values), in addition, both are written to the folder named 
##' "arsPathway".
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step (example 2)
##' #select.pathway.list from the pathSelect step
##' dataset.names <- c("hb","hs","ht","mb","ms","mt")
##' ARS_pathway <- multiARS_pathway(mcmc.merge.list,dataset.names,
##' select.pathway.list,B=100)
##' }

multiARS_pathway <- function(mcmc.merge.list,dataset.names,
                             select.pathway.list,
                             measure="Fmeasure",
                             parallel=F,cpu=2,B=50){
  
  names(mcmc.merge.list) <- dataset.names
  M <- length(mcmc.merge.list)
  P <- choose(M,2)
  select.pathways <- names(select.pathway.list)
  data_genes <- rownames(mcmc.merge.list[[1]])
  pathway.size <- sapply(select.pathway.list,function(x) {
                         length(intersect(data_genes,x))})
  K <- length(select.pathways)
  
  ARS <- ARSpvalue <- array(1,dim=c(K,M,M),dimnames=
                              list(select.pathways,dataset.names,dataset.names))

  for(k in 1:K){
    print(paste("pathway:",k,sep=""))
    path_genes <- select.pathway.list[[k]]
    for(i in 1:(M-1)){
      for(j in (i+1):M){
      
      dat1 <- mcmc.merge.list[[i]]
      dat2 <- mcmc.merge.list[[j]]
      deIndex1 <- attr(dat1,"DEindex") 
      deIndex2 <- attr(dat2,"DEindex")         
      names(deIndex1) <- rownames(dat1)[deIndex1]
      names(deIndex2) <- rownames(dat2)[deIndex2]
      data_genes <- rownames(dat1)
      
      genek <- intersect(path_genes,data_genes)
      dat1_k <- dat1[genek,]
      dat2_k <- dat2[genek,]
    
    if(length(intersect(names(deIndex1), genek))<=3 ){
      deIndex1_k <- 1:nrow(dat1_k)
    } else {
      deIndex1_k <- which(rownames(dat1_k)%in%intersect(names(deIndex1), genek))
    } 
    
    if(length(intersect(names(deIndex2), genek))<=3 ){
      deIndex2_k <- 1:nrow(dat2_k)
    } else {
      deIndex2_k <- which(rownames(dat2_k)%in%intersect(names(deIndex2), genek))
    } 
    
    ARS[k,i,j] <- ARS[k,j,i] <- ars(dat1_k,dat2_k,deIndex1_k,deIndex2_k,measure="Fmeasure")
  
      }
    }
}
  
  ARSperm <- array(0,dim=c(B,K,M,M))
  
  if(parallel==F){
    for(i in 1:(M-1)){
      for(j in (i+1):M){
        ARSperm[,,i,j] <- ARSperm[,,j,i] <- permPathway(dat1,dat2,select.pathways,pathway.size,B)
        print(paste("permuted part - pair: dataset ",i," and dataset ",j,sep=""))
      }
    }  
  } else if(parallel==T){
    for(i in 1:(M-1)){
      for(j in (i+1):M){
        ARSperm[,,i,j] <- ARSperm[,,j,i] <- permParPathway(dat1,dat2,select.pathways,pathway.size,B,cpu)
        print(paste("permuted part - pair: dataset ",i," and dataset ",j,sep=""))
      }
    }
  } 
  
  for(k in 1:K){
    for(i in 1:(M-1)){
      for(j in i:M){
        ARSperm_kij <- ARSperm[,k,i,j]
        ARSpvalue[k,i,j] <- ARSpvalue[k,j,i] <- arsPvalue(ARS[k,i,j],ARSperm_kij) 
      }
    }  
  }
  
  ARS.mat <- ARSpvalue.mat <- matrix(0,K,P)
  rownames(ARS.mat) <- rownames(ARSpvalue.mat) <- select.pathways
  
  combinations <- combn(dataset.names,m=2)
  colnames(ARS.mat) <- colnames(ARSpvalue.mat) <- apply(combinations,2,FUN=function(x) paste(x,collapse ="_"))
  
  ARS.mat <- data.frame(ARS.mat)
  ARSpvalue.mat <- data.frame(ARSpvalue.mat)
  
  for(p in 1:P){
    pairname <- colnames(ARS.mat)[p]
    name1 <- strsplit(pairname,split="_",fixed=T)[[1]][1]
    name2 <- strsplit(pairname,split="_",fixed=T)[[1]][2]
    ARS.mat[,pairname] <- ARS[,name1,name2]
    ARSpvalue.mat[,pairname] <- ARSpvalue[,name1,name2]
  }
  
  dir.path <- "arsPathway"
  if (!file.exists(dir.path)) dir.create(dir.path)
  write.csv(ARS.mat,file=paste(paste(dir.path,"ARS_pathway_",sep="/"),M,".csv",sep=""))
  write.csv(ARSpvalue.mat,file=paste(paste(dir.path,"ARSpvalue_pathway_",sep="/"),M,".csv",sep=""))
  out <- list(ARS.mat=ARS.mat,ARSpvalue.mat=ARSpvalue.mat)
  return(out)
}




