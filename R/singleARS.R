##' Resemblance analysis for single pair: global ARS and its permuted p-value
##' The \code{singleARS_global} is function to perform resemblance analysis 
##' for single pair, generating global ARS and its permuted p-value. 
##' @title Resemblance analysis for single pair: global ARS and its permuted 
##' p-value.
##' @param mcmc.merge.list: a list of merged MCMC output matrices.
##' @param measure: three types of ARS measures to be used: "youden",
##' "Fmeasure","geo.mean". Default is "Fmeasure". 
##' @param parallel: whether to perform parallel computing in permutation.
##' @param cpu: if parallel=T, how many cpus to be used.
##' @param B: number of permutations. 

##' @return global ARS values and its permuted p-value, in addition, the two 
##' values are written to the folder named "arsGlobal".
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step
##' ARS_global <- singleARS_global(mcmc.merge.list,B=100)
##' }

singleARS_global <- function(mcmc.merge.list,
                             measure="Fmeasure",parallel=F,cpu=2,B=50){
  
  dat1 <- mcmc.merge.list[[1]]
  dat2 <- mcmc.merge.list[[2]]
  deIndex1 <- attr(dat1,"DEindex") 
  deIndex2 <- attr(dat2,"DEindex") 
  
  if(parallel==F){
   ARS <- ars(dat1,dat2,deIndex1,deIndex2,measure=measure)
   ARSperm <- permGlobal(dat1,dat2,B=B)
   ARSpvalue <- arsPvalue(ARS,ARSperm)
   out <- c(ARS=ARS,ARSpvalue=ARSpvalue)
   dir.path <- "arsGlobal"
   if (!file.exists(dir.path)) dir.create(dir.path)
   write.csv(out,file=paste(dir.path,"ARS_global_single.csv",sep="/"))
  } else if(parallel==T){
    ARS <- ars(dat1,dat2,deIndex1,deIndex2,measure=measure)
    ARSperm <- permParGlobal(dat1,dat2,B=B,cpu=cpu)
    ARSpvalue <- arsPvalue(ARS,ARSperm)
    dir.path <- "arsGlobal"
    out <- c(ARS=ARS,ARSpvalue=ARSpvalue)
    dir.path <- "arsGlobal"
    if (!file.exists(dir.path)) dir.create(dir.path)
    write.csv(out,file=paste(dir.path,"ARS_global_single.csv",sep="/"))
  } 
  
  return(out)
}

##' Resemblance analysis for single pair: pathway specific ARS and their 
##' permuted p-value
##' The \code{singleARS_pathway} is function to perform resemblance analysis 
##' for single pair, generating pathway specific ARS and their permuted 
##' p-value. 
##' @title Resemblance analysis for single pair: pathway specific ARS and 
##' their permuted p-value.
##' @param mcmc.merge.list: a list of merged MCMC output matrices.
##' @param select.pathway.list: a list of selected pathways (containing gene 
##' components).
##' @param measure: three types of ARS measures to be used: "youden",
##' "Fmeasure","geo.mean". Default is "Fmeasure".
##' @param parallel: whether to perform parallel computing in permutation.
##' @param cpu: if parallel=T, how many cpus to be used.
##' @param B: number of permutations. 

##' @return a data frame of pathway specific ARS values and their permuted 
##' p-value (pathway on rows, 1st column being ARS value and 2nd column being
##' the p-values), in addition, the dataframe is written to the folder named 
##' "arsPathway".
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step
##' #select.pathway.list from the pathSelect step
##' ARS_pathway <- singleARS_pathway(mcmc.merge.list,select.pathway.list,B=100)
##' }

singleARS_pathway <- function(mcmc.merge.list,select.pathway.list,
                              measure="Fmeasure",
                              parallel=F,cpu=2,B=50){
  
  dat1 <- mcmc.merge.list[[1]]
  dat2 <- mcmc.merge.list[[2]]
  deIndex1 <- attr(dat1,"DEindex") 
  deIndex2 <- attr(dat2,"DEindex") 
  names(deIndex1) <- rownames(dat1)[deIndex1]
  names(deIndex2) <- rownames(dat2)[deIndex2]
  data_genes <- rownames(dat1)
  
  select.pathways <- names(select.pathway.list)
  pathway.size <- sapply(select.pathway.list,function(x) {
                            length(intersect(data_genes,x))})
  K <- length(select.pathways)
  ARS <- ARSpvalue <- rep(NA,K)
  for(k in 1:K){
    path_genes <- select.pathway.list[[k]]
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
    
    ARS[k] <- ars(dat1_k,dat2_k,deIndex1_k,deIndex2_k,measure="Fmeasure")
    print(k)
  }
  
  if(parallel==F){
    ARSperm <- permPathway(dat1,dat2,select.pathways,pathway.size,B)
  } else if(parallel==T){
    ARSperm <- permParPathway(dat1,dat2,select.pathways, pathway.size,B,cpu)
  } 
  
  for(k in 1:K){
    ARSperm_k <- ARSperm[,k]
    ARSpvalue[k] <- arsPvalue(ARS[k],ARSperm_k) 
  }
  
  out <- data.frame(ARS=ARS,ARSpvalue=ARSpvalue)
  rownames(out) <- select.pathways
  
  dir.path <- "arsPathway"
  if (!file.exists(dir.path)) dir.create(dir.path)
  write.csv(out,file=paste(dir.path,"ARS_pathway_single.csv",sep="/"))
  
  return(out)
}

ars <- function(dat1,dat2,deIndex1,deIndex2,
                measure=c("youden","Fmeasure","geo.mean")){
  ## same for both global and pathway
  if(measure=="youden"){
    ARS <- ARSY(dat1,dat2,deIndex1)
  } else if(measure=="Fmeasure"){
    ARS <- ARSF(dat1,dat2,deIndex1,deIndex2)
  } else if(measure=="geo.mean"){
    ARS <- ARSG(dat1,dat2,deIndex1)
  }
  return(ARS)
}

permGlobal <- function(delta1,delta2,B){
  G <- nrow(delta1)
  top <- 500
  out <- rep(NA,B)
  for(b in 1:B){
    delta1perm <- delta1[sample(1:nrow(delta1),G,replace = F),]
    delta2perm <- delta2[sample(1:nrow(delta2),G,replace = F),]
    permDE1 <- sample(1:G,top,replace = F)
    permDE2 <- sample(1:G,top,replace = F)
    out[b] <- round(ARSF(delta1perm,delta2perm,
                         permDE1,permDE2),digits=3)
  }  
  return(out)
}

permPathway <- function(delta1,delta2,
                        select.pathways, pathway.size,B){
  K <- length(select.pathways)
  G <- nrow(delta1)
  
  out <- matrix(NA,nrow=B,ncol=K)
  colnames(out) <- select.pathways
  
  for(b in 1:B){
    delta1perm <- delta1[sample(1:nrow(delta1),nrow(delta1),replace = F),]
    delta2perm <- delta2[sample(1:nrow(delta2),nrow(delta2),replace = F),]
    
    out[b,] <- sapply(1:K,function(k){
      size <- pathway.size[k]
      index <- sample(1:G,size,replace=F)
      d1.select <- delta1perm[index,]
      d2.select <- delta2perm[index,]
      if(length(index)<=5){
        d1.permDE <- d2.permDE <- 1:length(index)
      } else {
        d1.permDE <- sample(1:length(index),round(0.5*size),replace = F)
        d2.permDE <- sample(1:length(index),round(0.5*size),replace = F)
      }
      kout <- round(ARSF(d1.select,d2.select,
                         d1.permDE,d2.permDE),digits=3)
      return(kout)
    },simplify = T)  ## each x K matrix
  }  
  
  return(out)
}

permParGlobal <- function(delta1,delta2,B,cpu){
  each <- B/cpu
  G <- nrow(delta1)
  top <- 500
  
  sfInit(parallel=T,cpus=cpu,type="SOCK")
  
  parFun <- function(xxx) {
    
    delta1perm <- delta1[sample(1:nrow(delta1),G,replace = F),]
    delta2perm <- delta2[sample(1:nrow(delta2),G,replace = F),]
    permDE1 <- sample(1:G,top,replace = F)
    permDE2 <- sample(1:G,top,replace = F)
    
    out <- round(ARSF(delta1perm,delta2perm,
                      permDE1,permDE2),digits=3)
    return(out)
  }  
  
  sfExport(list=c("delta1","delta2","G","top"))
  
  result<-sfLapply(1:cpu, parFun) 
  sfStop()
  
  out <- unlist(result)
  return(out)
}

permParPathway <- function(delta1,delta2,
                           select.pathways, pathway.size,
                           B,cpu){
  each <- B/cpu
  K <- length(select.pathways)
  G <- nrow(delta1)
  
  sfInit(parallel=T,cpus=cpu,type="SOCK")
  
  parFun <- function(xxx) {
    
    delta1perm <- delta1[sample(1:nrow(delta1),nrow(delta1),replace = F),]
    delta2perm <- delta2[sample(1:nrow(delta2),nrow(delta2),replace = F),]
    
    out <- t(replicate(each,sapply(1:K,function(k){
      size <- pathway.size[k]
      index <- sample(1:G,size,replace=F)
      d1.select <- delta1perm[index,]
      d2.select <- delta2perm[index,]
      d1.permDE <- sample(1:length(index),round(0.2*size),replace = F)
      d2.permDE <- sample(1:length(index),round(0.2*size),replace = F)
      
      kout <- round(ARSF(d1.select,d2.select,
                         d1.permDE,d2.permDE),digits=3)
      return(kout)
    },simplify = T)))  ## each x K matrix
    
    colnames(out) <- select.pathways
    return(out)
  }  
  
  sfExport(list=c("delta1","delta2","select.pathways",
                  "pathway.size","G","K","each"))
  
  result<-sfLapply(1:cpu, parFun) 
  sfStop()
  
  out <- do.call(rbind, result)
  return(out)
}

arsPvalue <- function(ARS, ARSperm){
  ARSp <- (sum(ARSperm>=ARS) + 1)/(length(ARSperm)+1)
  return(ARSp)
}


