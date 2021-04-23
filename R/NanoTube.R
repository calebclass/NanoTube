#' NanoTube.
#' 
#' A package for NanoString nCounter gene expression data processing and
#' analysis.
#' 
#' @name NanoTube
#' @docType package
#' @import Biobase ggplot2 fgsea limma reshape
NULL

#' Example pathway database
#' 
#' A list object containing example gene sets from WikiPathways.
#' 
#' @name ExamplePathways
#' @docType data
#' @keywords datasets
#' @usage data(ExamplePathways)
#' @format A list object with 30 vectors of gene symbols, for 30 pathways
NULL

#' Example results from runLimmaAnalysis
#' 
#' Results of runLimmaAnalysis using the example data set GSE117751 (in extdata).
#' 
#' @name ExampleResults
#' @docType data
#' @keywords datasets
#' @usage data(ExampleResults)
#' @format An MArrayLM object from limma
NULL