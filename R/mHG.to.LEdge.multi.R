#' Process NanoString nCounter gene expression data.
#'
#' This function identifies leading edge genes from an mHG analysis result.
#'
#' @param mHG.res Result from mHG.genesets.multi
#' @param limmaResults t matrix from runLimmaAnalysis result
#' @param t.thresh t-statistic threshold to include gene in leading edge
#' @param q.thresh q-value threshold to include pathway in analysis
#' @param geneset.file Same geneset file used in mHG.genesets.multi


mHG.to.LEdge.multi <- function(mHG.res, limmaResults, t.thresh, q.thresh, geneset.file,
                               dict=NULL, colnames_length=NULL) {
#Remove dict, colnames_length!

  file.type <- substr(geneset.file, start = nchar(geneset.file)-2, stop = nchar(geneset.file))
  if (file.type == "rds"){
    genesets <- readRDS(geneset.file)
  } else if (file.type == "gmt"){
    genesets <- read.gmt(geneset.file)
  } else {
    stop("geneset.file must be in .rds (list object) or .gmt format")
  }

  le.results <- list()
  for (comp in names(mHG.res)){
    le.results[[comp]] <- list(positive = NULL,
                               negative = NULL)

    # Positive gene sets
    if (class(try(mHG.res[[comp]][[1]])) != "try-error") {
      if(any(mHG.res[[comp]][[1]]$p.adj < q.thresh & !is.na(mHG.res[[comp]][[1]]$p.adj))) {
        genesets.signif.pos <- genesets[rownames(mHG.res[[comp]][[1]])[mHG.res[[comp]][[1]]$p.adj < q.thresh & !is.na(mHG.res[[comp]][[1]]$p.adj)]]
        le.results[[comp]][[1]] <- matrix(0, nrow = sum(limmaResults$t[,comp] > t.thresh), ncol = length(genesets.signif.pos),
                                          dimnames = list(rownames(limmaResults$t)[limmaResults$t[,comp] > t.thresh], names(genesets.signif.pos)))
        for (set in names(genesets.signif.pos)){
          le.results[[comp]][[1]][rownames(le.results[[comp]][[1]]) %in% genesets.signif.pos[[set]],set] <- 1
        }
        if (!is.null(dict)) rownames(le.results[[comp]][[1]]) <- dict$Name[match(rownames(le.results[[comp]][[1]]), dict$genes.entrez)]
        if (!is.null(colnames_length)) colnames(le.results[[comp]][[1]]) <- substr(colnames(le.results[[comp]][[1]]), 1, colnames_length)
        if (ncol(le.results[[comp]][[1]]) > 1) le.results[[comp]][[1]] <- le.results[[comp]][[1]][rowSums(le.results[[comp]][[1]]) > 0,]
      }
    }

    # Negative gene sets
    if (class(try(mHG.res[[comp]][[2]])) != "try-error") {
      if(any(mHG.res[[comp]][[2]]$p.adj < q.thresh & !is.na(mHG.res[[comp]][[2]]$p.adj))) {
        genesets.signif.neg <- genesets[rownames(mHG.res[[comp]][[2]])[mHG.res[[comp]][[2]]$p.adj < q.thresh & !is.na(mHG.res[[comp]][[2]]$p.adj)]]
        le.results[[comp]][[2]] <- matrix(0, nrow = sum(limmaResults$t[,comp] < -t.thresh), ncol = length(genesets.signif.neg),
                                          dimnames = list(rownames(limmaResults$t)[limmaResults$t[,comp] < -t.thresh], names(genesets.signif.neg)))
        for (set in names(genesets.signif.neg)){
          le.results[[comp]][[2]][rownames(le.results[[comp]][[2]]) %in% genesets.signif.neg[[set]],set] <- 1
        }
        if (!is.null(dict)) rownames(le.results[[comp]][[2]]) <- dict$Name[match(rownames(le.results[[comp]][[2]]), dict$genes.entrez)]
        if (!is.null(colnames_length)) colnames(le.results[[comp]][[2]]) <- substr(colnames(le.results[[comp]][[2]]), 1, colnames_length)
        if (ncol(le.results[[comp]][[2]]) > 1) le.results[[comp]][[2]] <- le.results[[comp]][[2]][rowSums(le.results[[comp]][[2]]) > 0,]
      }
    }
  }

  return(le.results)
}
