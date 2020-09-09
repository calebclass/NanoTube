#' Run all gene set analyses
#'
#' Runs minimum hypergeometric tests using the mHG library. This conducts mHG.genesets on every
#' possible comparison from the Limma analysis, in both the enriched (positive t-statistics) and
#' decreased (negative t-statistics) directions.
#'
#' @param limmaResults Result from runLimmaAnalysis.
#' @param geneset.file Gene set file name, in .rds (list), .gmt, or .tab format. Gene names must be
#' in the same form as in the ranked.list.
#' @param sourceDB Source database to include, only if using a .tab-format 
#' geneset.file from CPDB.
#' @param t_thresh The threshold t-statistic for considering a gene positively or negatively.
#' For "negative" analysis, this threshold is set to -t_thresh.
#' @param numReq_overall The number of genes that must be differentially expressed above t_thresh
#' to conduct mHG analysis. If this is not satisfied, the comparison will be skipped.
#' @param numReq_inSet Number of genes required to conduct analysis on a given gene set (default = 1).
#' If fewer than this number of genes from the ranked list are included in a gene set, that
#' gene set will be skipped for this analysis.
#' @return A list containing data frames with the mHG results for each comparison. Within each
#' comparison, both a "positive" and "negative" results data frame are provided, with mHG results
#' for genes most enriched in the comparison group vs. the base group (positive), and most enriched
#' in the base group vs. the comparison group (negative).


mHG.genesets.multi <- function(limmaResults, geneset.file, sourceDB = NULL,
                               t_thresh = 2, numReq_overall = 10, numReq_inSet = 5,
                               skipFirst = TRUE) {
  # limmaResults -- Results from runLimmaAnalysis
  # geneset.file -- geneset file name
  # t_thresh -- t-stat minimum for inclusion in mHG test
  # numReq_overall -- minimum number of genes above t_thresh to conduct mHG analysis on a factor
  # numReq_inSet -- minimum number of genes from geneset in codeset to analyze a given geneset

  mHG.res <- list()

  for (i in ifelse(skipFirst, yes=2, no=1):ncol(limmaResults$t)) {
    testDat <- limmaResults$t[,i]
    testDat <- testDat[order(testDat, decreasing = TRUE)]
    mHG.res[[colnames(limmaResults$t)[i]]] <- list()
    if (sum(testDat > t_thresh) >= numReq_overall) {
      mHG.res[[colnames(limmaResults$t)[i]]]$positive <-
        mHG.genesets(testDat,
                     geneset.file, sourceDB,
                     numReq = numReq_inSet, n_max = sum(testDat > t_thresh))
    } else {
      mHG.res[[colnames(limmaResults$t)[i]]]$positive <- NULL
    }
    testDat2 <- rev(-testDat)
    if (sum(testDat2 > t_thresh) >= numReq_overall) {
      mHG.res[[colnames(limmaResults$t)[i]]]$negative <-
        mHG.genesets(testDat2,
                     geneset.file, sourceDB,
                     numReq = numReq_inSet, n_max = sum(testDat2 > t_thresh))
    } else {
      mHG.res[[colnames(limmaResults$t)[i]]]$negative <- NULL
    }
  }

  return(mHG.res)
}
