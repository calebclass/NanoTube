#' Run gene set enrichment analysis using DE results.
#' 
#' Use the fgsea library to run gene set enrichment analysis from the Limma
#' analysis results. Genes will be ranked by their log2 fold changes
#' or t-statistics.
#' 
#' @param limmaResults Result from runLimmaAnalysis.
#' @param gene.sets Gene set file name, in .rds (list), .gmt, or .tab format;
#' or a list object containing the gene sets. Gene names must be
#' in the same form as in the ranked.list.
#' @param sourceDB Source database to include, only if using a .tab-format 
#' geneset.file from CPDB.
#' @param min.set Number of genes required to conduct analysis on a given gene 
#' set (default = 1). If fewer than this number of genes from limmaResults are 
#' included in a gene set, that gene set will be skipped for this analysis.
#' @param rank.by Rank genes by log2 fold changes ('coefficients', default) or 
#' t-statistics ('t').
#' @param skip.first Logical: Skip the first factor for gene set analysis?
#' Frequently the first factor is the 'Intercept', which is generally 
#' uninteresting for GSEA (default TRUE).
#' @return A list containing data frames with the fgsea results for each 
#' comparison. 

limmaToFGSEA <- function(limmaResults, gene.sets, sourceDB = NULL,
                         min.set = 1, rank.by = c('coefficients', 't'),
                         skip.first = TRUE) {

  if (class(gene.sets) == "list") {
    gene.set.list <- gene.sets
  } else {
    file.type <- substr(gene.sets, start = nchar(gene.sets)-2, stop = nchar(gene.sets))
    if (file.type == "rds") {
      gene.set.list <- readRDS(gene.sets)
    } else if (file.type == "gmt") {
      gene.set.list <- read.gmt(gene.sets)
    } else if (file.type == "tab") {
      gene.set.list <- read_cpdb_tab(gene.sets, sourceDB)
    } else {
      stop("gene.sets must be in .rds (list object), .gmt, or .tab format")
    }
  }
  
  rank.by <- rank.by[1]
  if (!(rank.by) %in% c('coefficients', 't')) {
    stop("rank.by must be either 'coefficeints' or 't'.")
  }
  
  if (skip.first) {
    analysis.factors <- colnames(limmaResults[[rank.by]])[-1]
  } else {
    analysis.factors <- colnames(limmaResults[[rank.by]])
  }
  
  fgsea.res <- list()
  for (i in analysis.factors) {
    fgsea.res[[i]] <- fgsea::fgseaMultilevel(pathways = gene.set.list,
                                             stats = limmaResults[[rank.by]][,i],
                                             minSize = min.set)
  }
  
  return(fgsea.res)
  
}
  