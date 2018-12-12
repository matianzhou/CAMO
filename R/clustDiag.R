##' Perform consensus clustering diagnostics to determine optimal K for pathway 
##' clustering 
##' The \code{clustDiag} is function to perform consensus clustering diagnostics for 
##' pathway clustering. It will generate consensus CDF and delta area plot and saved to 
##' the folder named "clustDiag". These plots will help you determine the optimal number
##' of clusters K for pathway clustering. Run this first before \code{multiOutput} when 
##' "clustPathway" output is chosen. 
##' @title Perform consensus clustering diagnostics to determine optimal K for pathway 
##' clustering. 
##' @param ARS_pathway: a list of two data frames: pathway specific ARS values and 
##' their permuted p-value (pathway on rows, column being ARS value or the p-values).
##' 
##' @return stored output in the folder named "clustDiag".
##' @export
##' @examples
##' \dontrun{
##' #ARS_pathway from the multiARS step 
##' results <- clustDiag(ARS_pathway)
##' }


clustDiag <- function(ARS_pathway){
  ### Clustering diagnosis
  orig.path <- getwd()
  ARSpvalue.mat <- ARS_pathway[["ARSpvalue.mat"]]
  
  dir.path <- "clustDiag"
  if (!file.exists(dir.path)) dir.create(dir.path)
  setwd(paste(orig.path,"/",dir.path,sep=""))
  
  #consensus clustering diagnostics
  results <- clustPathway(ARSpvalue.mat)
  setwd(orig.path)
  return(results)
}