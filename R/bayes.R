##' Run Bayesian analysis for individual pData 
##' The \code{bayes} is function to run Bayesian analysis for individual pdata
##' @title Run Bayesian analysis for individual pdata
##' @param pData: individual pData. The data matrix has to consist of 
##' two columns, first being the p-value, second being the log fold change.
##' @param seed: seed number. 

##' @return an MCMC output matrix of signed DE indicator (input for
##' resemblance analysis) 
##' @export
##' @examples
##' \dontrun{
##' data(hb)
##' summaryDE <- indDE(data=data,group=group,data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' pData <- summaryDE[,c(3,1)]
##' MCMCout <- bayes(pData, seed=12345)
##' }


bayes <- function(pData, seed=12345){

  set.seed(seed)
  
## deSelect part
  
  DEindex <- deSelect(pData)
  
## bayesP part   
  p <- pData[,1]
  lfc <- pData[,2]
  G <- nrow(pData)
  z <- PtoZ(p,lfc)
  iteration <- 5000
  burnin <- 3000
  thin <- 10
  names(z) <- rownames(pData)
  prop <- SelectGamma(p)
  if(prop <= 0.3) {
    gamma <- G*0.3
  } else{
    gamma <- G*prop
  }
  MCMCout <- MCMC(z, iteration, gamma) 
  signdelta <- MCMCout$Y[,-c(1:burnin)]
  
  signdelta.sub <- signdelta[,seq(1,ncol(signdelta),by=thin)]
  rownames(signdelta.sub) <- names(z)
  
  attr(signdelta.sub,"DEindex") <- DEindex
  
  return(signdelta.sub) #the full sign delta for a dataset (subsample 500)
}  


deSelect <- function(pData, q.cut=0.3,topDE.number = 1000){
  pvalue <- pData[,1]
  qvalue <- p.adjust(pvalue,method="BH")
  if (is.null(q.cut) || sum(qvalue < q.cut) <=topDE.number ) { 
    DEindex <- which(rank(pvalue,ties.method="first") %in% 1:topDE.number)
    names(DEindex) <- names(pvalue)[DEindex]
  } else {
    DEindex <- which(qvalue < q.cut)
    names(DEindex) <- names(pvalue)[DEindex]
  }
  return(DEindex)
}