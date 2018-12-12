##########################
###### Check functions ###
##########################

check.rawData <- function(data){
  G <- nrow(data)
  N <- ncol(data)
  if(!is.numeric(data)) {
    stop("expression data not numeric")
  }
  if(is.null(row.names(data))) {
    stop("gene symbol missing")
  }
  if(N <= 3) {
    stop("too few samples")
  }
  if(sum(duplicated(row.names(data)))>0){
  	stop("duplicate gene symbols")
  } 
}

check.groupData <- function(group){
  l <- nlevels(group)
  if( l != 2 ) {
    stop("not a two-class comparison")
  }
}

check.pData <- function(pData){
  G <- nrow(pData)
  if(!is.numeric(pData)) {
    stop("pvalue data not numeric")
  }
  if(is.null(row.names(pData))) {
    stop("gene symbol missing")
  }
  if(ncol(pData) <2) {
    stop("missing either p-value or effect size")
  }
}


check.compatibility <- function(data, group, case.label, ctrl.label){
  G <- nrow(data)
  N1 <- ncol(data)
  N2 <- length(group)
  if(N1 != N2) {
    stop("expression data and class label have unmatched sample size")
  }
  if(!all(group %in% c(case.label,ctrl.label))){
    stop("including class labels other than the case and control")
  }
  if(sum(group==case.label) <= 1 ||sum(group==ctrl.label) <= 1){
    stop("not enough samples in either case or control group")
  }
}

##########################
###### BayesP part ###
##########################

PtoZ <- function(p2, lfc) {
  sgn <- sign(lfc)
#  z <- ifelse(sgn>0, qnorm(1-p2/2,lower.tail=T), qnorm(p2/2,lower.tail = T))
  z <- ifelse(sgn>0, qnorm(p2/2,lower.tail = F), qnorm(p2/2,lower.tail = T))
  return(z)
}

SelectGamma <- function(p){
  ## Gamma is DE proportion  = 1-pi0
  m <- length(p)
  lambda <- seq(0,0.95,by=0.01)
  pi0 <- sapply(lambda, function(x) sum(p>x)/(m*(1-x))  )
  
  # fit a natural cubic spline
  library(splines)
  dat <- data.frame(pi0=pi0, lambda=lambda)
  lfit <- lm(pi0 ~ ns(lambda, df = 3), data=dat)
  pi0hat <- predict(lfit, data.frame(lambda=1))
  gamma <- 1 - pi0hat
  return(gamma)
} 


##########################
###### ARS functions ###
##########################

ARSY <- function(d1,d2,index1){
  Sens <- calcSens(d1[index1,],d2[index1,])
  Spec <- calcSpec(d1[-index1,],d2[-index1,])
  ESens <- calcESens(d1,d2)
  ESpec <- calcESpec(d1,d2)
  RSY <- Sens + Spec - 1
  ERSY <- ESens + ESpec - 1
  ARS <- (RSY - ERSY)/(1-ERSY)
  if(is.nan(ARS)) {
    ARS <- 0
  }
  return(ARS)
}


ARSF <- function(d1,d2,index1,index2){
  Sens <- calcSens(d1[index1,],d2[index1,])
  Prec <- calcPrec(d1[index2,],d2[index2,])
  ESens <- calcESens(d1,d2)
  EPrec <- calcEPrec(d1,d2)
  RSF <- (2*Sens*Prec)/(Sens+Prec)
  ERSF <- (2*ESens*EPrec)/(ESens+EPrec)
  ARS<- (RSF - ERSF)/(1-ERSF)
  if(is.nan(ARS)) {
    ARS <- 0
  }
  return(ARS)
}

ARSG <- function(d1,d2,index1){
  Sens <- calcSens(d1[index1,],d2[index1,])
  Spec <- calcSpec(d1[-index1,],d2[-index1,])
  ESens <- calcESens(d1,d2)
  EPrec <- calcEPrec(d1,d2)
  ESpec <- calcESpec(d1,d2)
  RSG<- sqrt(Sens*Spec)
  ERSG <- sqrt(ESens*ESpec)
  ARS <- (RSG- ERSG)/(1-ERSG)
  if(is.nan(ARS)) {
    ARS <- 0
  }
  return(ARS)
}

ARSYperm <- function(d1,d2){
  Sens <- calcSens(d1,d2)
  Spec <- calcSpec(d1,d2)
  ESens <- calcESens(d1,d2)
  ESpec <- calcESpec(d1,d2)
  RSY <- Sens + Spec - 1
  ERSY <- ESens + ESpec - 1
  ARS <- (RSY - ERSY)/(1-ERSY)
  if(is.nan(ARS)) {
    ARS <- 0
  }
  return(ARS)
}


ARSFperm <- function(d1,d2){
  Sens <- calcSens(d1,d2)
  Prec <- calcPrec(d1,d2)
  ESens <- calcESens(d1,d2)
  EPrec <- calcEPrec(d1,d2)
  RSF <- (2*Sens*Prec)/(Sens+Prec)
  ERSF <- (2*ESens*EPrec)/(ESens+EPrec)
  ARS<- (RSF - ERSF)/(1-ERSF)
  if(is.nan(ARS)) {
    ARS <- 0
  }
  return(ARS)
}

ARSGperm <- function(d1,d2){
  Sens <- calcSens(d1,d2)
  Spec <- calcSpec(d1,d2)
  ESens <- calcESens(d1,d2)
  EPrec <- calcEPrec(d1,d2)
  ESpec <- calcESpec(d1,d2)
  RSG<- sqrt(Sens*Spec)
  ERSG <- sqrt(ESens*ESpec)
  ARS <- (RSG- ERSG)/(1-ERSG)
  if(is.nan(ARS)) {
    ARS <- 0
  }
  return(ARS)
}

ARStransform <- function(ARS, theta=7) {
  trun.ARS <-ifelse(ARS<0,0,ARS)
  trsf.ARS <- theta*exp(-theta*trun.ARS)
  return(trsf.ARS)
}

##########################
##Pathway enrich analysis###
##########################

gsa.fisher <- function(x, background, pathway) {
  ####x is the list of query genes
  ####backgroud is a list of background genes that query genes from 
  ####pathway is a list of different pathway genes
  count_table<-matrix(0,2,2)
  x<-toupper(x)
  background<-toupper(background)
  index<-which(toupper(background) %in% toupper(x)==FALSE)
  background_non_gene_list<-background[index]
  x<-toupper(x)
  pathway<-lapply(pathway,function(x) intersect(toupper(background),toupper(x)))
  get.fisher <- function(path) {
    res <- NA
    ####in the gene list and in the pathway
    count_table[1,1]<-sum(x %in% path)
    #count_table[1,1]<-sum(is.na(charmatch(x,path))==0)
    ####in the gene list but not in the pathway
    count_table[1,2]<-length(x)-count_table[1,1]
    ####not in the gene list but in the pathway
    count_table[2,1]<-sum(background_non_gene_list%in% path)
    ####not in the gene list and not in the pathway
    count_table[2,2]<-length(background_non_gene_list)-count_table[2,1]       
    matched_gene<-x[x %in% path]
    match_num<-length(matched_gene)
    overlap_info<-array(0,dim=4)
    names(overlap_info)<-c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
    overlap_info[1]=count_table[1,1]
    overlap_info[2]=count_table[1,2]
    overlap_info[3]=count_table[2,1]
    overlap_info[4]=count_table[2,2]
    if(length(count_table)==4){
      res <- fisher.test(count_table, alternative="greater")$p}
    return(list(p_value=res,match_gene=matched_gene,match_num=match_num,
                fisher_table=overlap_info))
  }
  p_val<-rep(0,length(pathway))
  
  match_gene_list<-list(length(pathway))
  match_gene <- array(0,dim=length(pathway))
  num1<-array(0,dim=length(pathway))
  num2<-matrix(0,nrow=length(pathway),ncol=4)
  colnames(num2)<-c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
  for(i in 1:length(pathway)){
    result<-get.fisher(pathway[[i]])
    p_val[i]<-result$p_value
    match_gene_list[[i]]<-result$match_gene
    match_gene[i]<-paste(match_gene_list[[i]],collapse="/")
    num1[i]<-result$match_num
    num2[i,]<-result$fisher_table
  }
  names(p_val) <- names(pathway)
  q_val <- p.adjust(p_val, "BH")
  
  summary<-data.frame(pvalue=p_val,
                      qvalue=q_val,
                      DE_in_Set=num2[,1],
                      DE_not_in_Set=num2[,2],
                      NonDE_in_Set=num2[,3],
                      NonDE_not_in_Set=num2[,4])
  
  a<-format(summary,digits=3)    
  
  return(a)
}

fisher <- function(x){
  n <- length(x)
  y <- -2*log(x)
  Tf <- sum(y)
  return(1-pchisq(Tf,2*n))
}


##########################
###### SA functions ###
##########################

## energy function 

E_tot <- function(delta.mat,a,delta.est){
  ## vector "a" of length n
  ## vector "delta_est" of length K+1: start from theta_0, then ordered from k=1 to K
  n <- nrow(delta.mat)
  K <- length(unique(a))
  #theta_0 <- delta.est[1]
  theta_0 <- 0
  E <- sum(sapply(1:n, function(x) {
    sum(sapply(1:n, function(y){
      if(a[x]==a[y]){
        (delta.mat[x,y] - delta.est[as.character(a[x])])^2
      } else{
        (delta.mat[x,y] - theta_0)^2
      }
    },simplify=T))    
  }, simplify=T))

  return(E)
}

## Estimate of delta.est (the means)

Est_mean <- function(delta.mat,a){
  # the first element is always the off-diagonal parts
  n <- nrow(delta.mat)
  K <- length(unique(a))
  total <- rep(0,K+1)
  size <- rep(0,K+1)
  deltamean <- rep(0,K+1)
  names(deltamean) <- c(0,sort(unique(a)))
  for(k in 1:K){
    a_k <- sort(unique(a))[k]
    total[1+k] <- sum(delta.mat[a==a_k,a==a_k])
    size[1+k] <- sum(a==a_k)^2
    deltamean[1+k] <- total[1+k]/size[1+k]
  }
  deltamean[1] <- (sum(delta.mat) - sum(total[-1]))/(n*n - sum(size[-1]))
  return(deltamean)
}

## Trial = split or relocate

Split <- function(a) {
  n <- length(a)
  ua <- unique(a)
  if(length(ua)==n) {
    return(a)
  } else{
    #ua.pick <- sample(x=ua,size=1)
    #a[names(sample(x=which(a==ua.pick),size=1))] <- max(ua)+1
    a.pick <- sample(x=a,size=1)
    pick.ind <- sample(x=which(a==a.pick),size=1)
    a[pick.ind] <- max(a)+1
    return(a)
  }  
}

Relocate <- function(a){
  n <- length(a)
  ua <- unique(a)
  if(length(ua)==1) {
    return(a)
  } else{
    pick.ind <- sample(x=1:n,size=1)
    #a.pick <- a[pick.ind]
    #a[pick.ind] <- sample(x=a[-which(a==a.pick)],size=1)
    a[pick.ind] <- sample(x=a[-pick.ind],size=1)
    return(a)
  }
}


##########################
###### Scatterness ###
##########################


scatter <- function(dat,cluster.assign) {
  ## dat: K pathways on rows and pairwise on columns
  C <- length(unique(cluster.assign))
  w <- 3 ## c(3,9)
  alpha <- 0.01 ## c(0.01,0.05,0.1,0.2)
  n <- nrow(dat)
  p <- ncol(dat)
  chi2 <- qchisq(p = alpha,df=p,ncp=0)
  yc <- uc <- matrix(NA,nrow=n,ncol=p)
  muc <- matrix(NA,nrow=C,ncol=p)
for(c in 1:C){
  datc <- dat[cluster.assign==c,]
  n <- nrow(datc)
  muc[c,] <- apply(datc,2,mean)
  yc[cluster.assign==c,] <- apply(datc,1,function(x) x - muc[c,])
}

s <- median(apply(yc,1,function(x) median(abs(x-median(x)))))
uc <- yc/(w*s)
index <- which(abs(uc)<1)
sbw <- sqrt(n*p)*sqrt(sum(c(yc)[index]^2*(1-c(uc)[index]^2)^4))/
  abs(sum( (1-c(uc)[index])*(1-5*c(uc)[index]^2)))
rK <- sbw*sqrt(chi2) ## common radius

logical.scatter <- rep(NA,nrow(yc))

for(i in 1:nrow(dat)){
  dati <- dat[i,]
  dev <- sapply(1:C, function(c) sqrt(sum(dati-muc[c,])^2))
  logical.scatter[i] <- ifelse(sum(dev>=rK)==C,T,F)
  print(i)
}

scatter.index <- which(logical.scatter==T)
return(scatter.index)

}

##########################
###### Text mining #######
##########################

TextMine <- function(hashtb, pathways, pathway, result){
  k <- length(unique(result))
  cat("Performing Text Mining Analysis...\n")
  hashtb = hashtb[hashtb [,2]%in%which(pathways %in% pathway),]
  tmk = list()
  nperm = 1000
  if (nrow(hashtb) == 0){
    for (i in 1:(k-1)){
      tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)
    }
  }
  else{
    for (i in 1:k){
      e = result[result == i]
      e = which(pathways %in% names(e))
      hashcl = hashtb[hashtb [,2]%in%e,]
      hashcl = hashcl[duplicated(hashcl[,1]) | duplicated(hashcl[,1], fromLast=TRUE),]
      if (nrow(hashcl) != 0){
        hashf = hashcl
        hashf[,1] = 1
        hashf = aggregate(hashf[,-2] ~ rownames(hashf),data=hashf, FUN=sum)
        rownames(hashf) = hashf[,1]
        hashf = hashf[,-1]
        colnames(hashf) = c("count","sum")
        hashap = hashcl
        hashap[] = 0 
        mperm = matrix(nrow = nrow(hashf),ncol = nperm)
        for (j in 1:nperm){
          subtb = hashtb[hashtb [,2]%in%sample(1:length(pathways),length(e)),]
          subtb = rbind(subtb,hashap)
          subtb = subtb[rownames(subtb) %in% rownames(hashap),]
          subtb[,1] = 1
          subtb = aggregate(subtb[,-2] ~ rownames(subtb),data=subtb, FUN=sum)
          rownames(subtb) = subtb[,1]
          subtb = subtb[,-1]
          colnames(subtb) = c("count","sum")
          subtb = subtb[rownames(subtb) %in% rownames(hashf),]
          mperm[,j] = subtb[,2]
        }
        hashf[,"p-value"] = apply(cbind(hashf[,2],mperm),1,
                                  function(x)((nperm + 2)-rank(x)[1])/(nperm + 1))
        hashf[,"q-vlaue"] = p.adjust(hashf[,"p-value"],method = "BH")
        tmk[[i]] = hashf[order(hashf[,3],-hashf[,2]),]
      }
      else {tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)}
    }
  }
  return(tmk)
  
} # End of Text Mining


writeTextOut <- function(tm_filtered,k,pathway.summary) {
  cat("Cluster 1\n", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
  cat("Key words,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
  write.table(t(rownames(tm_filtered[[1]])[1:15]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, 
              append = T, row.names=F,col.names=F,na="")
  cat("q_value,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
  write.table(t(tm_filtered[[1]][1:15,4]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, 
              append = T, row.names=F,col.names=F,na="")
  cat("count,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
  write.table(t(tm_filtered[[1]][1:15,1]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, 
              append = T, row.names=F,col.names=F,na="")
  write.table(pathway.summary[[1]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, 
              append = T, row.names=F,col.names=F)
for (i in 2:k){ 
    cat(paste("\nCluster ", i, "\n", sep = ""), file = paste("Clustering_Summary_K",k,".csv",sep=""), append = T)
    cat("Key words,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(rownames(tm_filtered[[i]])[1:15]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, 
                append = T, row.names=F,col.names=F,na="")
    cat("q_value,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(tm_filtered[[i]][1:15,4]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, 
                append = T, row.names=F,col.names=F,na="")
    cat("count,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(tm_filtered[[i]][1:15,1]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, 
                append = T, row.names=F,col.names=F,na="")
    write.table(pathway.summary[[i]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, 
                append = T, row.names=F,col.names=F)
    
  }
}


