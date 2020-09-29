#' Conduct leading edge analysis for gene set analysis results
#'
#' This function identifies leading edge genes from an mHG analysis result. This is mostly
#' redundant with mHG.to.LEdge.multi.
#'
#' @param mHG.res Result from mHG.genesets
#' @param t.stats t-stats from analysis
#' @param t.thresh t-statistic threshold to include gene in leading edge
#' @param q.thresh q-value threshold to include pathway in analysis
#' @param geneset.file Same geneset file used in mHG.genesets.multi


mHG.to.LEdge <- function(mHG.res, t.stats, t.thresh, q.thresh, geneset.file, dict=NULL,
                         colnames_length = NULL) {

  file.type <- substr(geneset.file, start = nchar(geneset.file)-2, stop = nchar(geneset.file))
  if (file.type == "rds"){
    genesets <- readRDS(geneset.file)
  } else if (file.type == "gmt"){
    genesets <- read.gmt(geneset.file)
  } else {
    stop("geneset.file must be in .rds (list object) or .gmt format")
  }

  le.results <- list()

  # Positive gene sets
  if (class(try(mHG.res[[1]])) != "try-error") {
    if(any(mHG.res[[1]]$p.adj < q.thresh & !is.na(mHG.res[[1]]$p.adj))) {
      genesets.signif.pos <- genesets[rownames(mHG.res[[1]])[mHG.res[[1]]$p.adj < q.thresh & !is.na(mHG.res[[1]]$p.adj)]]
      le.results[[1]] <- matrix(0, nrow = sum(t.stats > t.thresh), ncol = length(genesets.signif.pos),
                                        dimnames = list(names(t.stats)[t.stats > t.thresh], names(genesets.signif.pos)))
      for (set in names(genesets.signif.pos)){
        le.results[[1]][rownames(le.results[[1]]) %in% genesets.signif.pos[[set]],set] <- 1
      }
      if (!is.null(dict)) rownames(le.results[[1]]) <- dict$Name[match(rownames(le.results[[1]]), dict$genes.entrez)]
      if (!is.null(colnames_length)) colnames(le.results[[1]]) <- substr(colnames(le.results[[1]]), 1, colnames_length)
      if (ncol(le.results[[1]]) > 1) le.results[[1]] <- le.results[[1]][rowSums(le.results[[1]]) > 0,]
    }}

  # Negative gene sets
  if (class(try(mHG.res[[2]])) != "try-error") {
    if(any(mHG.res[[2]]$p.adj < q.thresh & !is.na(mHG.res[[2]]$p.adj))) {
      genesets.signif.neg <- genesets[rownames(mHG.res[[2]])[mHG.res[[2]]$p.adj < q.thresh & !is.na(mHG.res[[2]]$p.adj)]]
      le.results[[2]] <- matrix(0, nrow = sum(t.stats < -t.thresh), ncol = length(genesets.signif.neg),
                                        dimnames = list(names(t.stats)[t.stats < -t.thresh], names(genesets.signif.neg)))
      for (set in names(genesets.signif.neg)){
        le.results[[2]][rownames(le.results[[2]]) %in% genesets.signif.neg[[set]],set] <- 1
      }
      if (!is.null(dict)) rownames(le.results[[2]]) <- dict$Name[match(rownames(le.results[[2]]), dict$genes.entrez)]
      if (!is.null(colnames_length)) colnames(le.results[[2]]) <- substr(colnames(le.results[[2]]), 1, colnames_length)
      if (ncol(le.results[[2]]) > 1) le.results[[2]] <- le.results[[2]][rowSums(le.results[[2]]) > 0,]
    }}

  return(le.results)

}
