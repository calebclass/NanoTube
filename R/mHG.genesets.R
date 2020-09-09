#' Run single gene set analysis.
#'
#' Runs a minimum hypergeometric test on the ranked.list for the genes in every
#' gene set of the provided file, using the mHG library. The ranked.list is analyzed
#' from 1 to n_max. If analyses in both directions (beginning and end of list, or positively-
#' and negatively-enriched genes) are desired, mHG.genesets.multi can be used.
#'
#' @param ranked.list A named vector, ranked in order of differential expression.
#' @param geneset.file Gene set file name, in .rds (list), .gmt, or .tab format. Gene names must be
#' in the same form as in the ranked.list.
#' @param sourceDB Source database to include, only if using a .tab-format 
#' geneset.file from CPDB.
#' @param numReq Number of genes required to conduct analysis on a given gene set (default = 1).
#' If fewer than this number of genes from the ranked list are included in a gene set, that
#' gene set will be skipped for this analysis.
#' @param n_max Maximum index for mHG test. By default, will allow the mHG test to include any
#' amount of genes from the list. This parameter can be adjusted to provide a maximum (e.g. the
#' first half of the list) length to be considered.
#' @return A data frame containing the mHG results for each gene set.

mHG.genesets <- function(ranked.list, geneset.file, sourceDB = NULL,
                         numReq = 1, n_max = length(ranked.list)) {
  library(mHG)

  file.type <- substr(geneset.file, start = nchar(geneset.file)-2, stop = nchar(geneset.file))
  if (file.type == "rds") {
    genesets <- readRDS(geneset.file)
  } else if (file.type == "gmt") {
    genesets <- read.gmt(geneset.file)
  } else if (file.type == "tab") {
    genesets <- read_cpdb_tab(geneset.file, sourceDB)
  } else {
    stop("geneset.file must be in .rds (list object) or .gmt format")
  }

  mHG.results <- as.data.frame(matrix(NA, nrow = length(genesets), ncol = 5,
                                         dimnames = list(names(genesets), c("numSet", "numCodeset", "mHG.index", "mHG.sum", "p.val"))))
  for (i in 1:length(genesets)){
    if (sum(names(ranked.list) %in% genesets[[i]]) >= numReq) {
      testDat <- ifelse(names(ranked.list) %in% genesets[[i]], yes = 1, no = 0)
      test <- mHG.test(testDat, n_max)
      mHG.results[i,] <- c(length(genesets[[i]]), sum(testDat),
                           test$n, test$b, test$p.value)
    }
  }

  mHG.results$p.adj <- p.adjust(mHG.results$p.val, method = "BH")

  return(mHG.results)
}
