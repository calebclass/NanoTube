#' Make master table of all mHG results
#'
#' This function clusters mHG results by leading edge similarity, and then combines to a
#' data frame or text file.
#'
#' @param genesetResults Results from pathway analysis using mHG.geneset.multi.
#' @param leadingEdge Results from mHG.to.LEdge.multi.
#' @param join.threshold The threshold distance to join gene sets. Gene sets with a distance
#' below this value will be joined to a single "cluster."
#' @param filename File name for the output text file. If NULL (default), data will be returned
#' as an R data frame.

make.mHG.master.table <- function(genesetResults, leadingEdge,
                                  join.threshold = 0.5,
                                  filename = NULL) {

  geneset.clustered <- list()

  for (i in names(genesetResults)) {
    for (j in names(genesetResults[[i]])) {
      nm <- paste(i, j, sep="-")

      if (is.null(leadingEdge[[i]][[j]])) {
        temp <-  cbind(rownames(genesetResults[[i]][[j]]), genesetResults[[i]][[j]][,-1])
        temp$Cluster <- NA
        temp$Cluster.Max <- ""
        colnames(temp)[c(1:4,8)] <- c("gene set","genes", "index", "count", "best")
        geneset.clustered[[nm]] <- temp[,c(7,8,1,5,6,2:4)]
      } else {
        geneset.clustered[[nm]] <- groupByDistance(genesetResults[[i]][[j]], leadingEdge[[i]][[j]],
                                                       join.threshold, returns="all")
      }

      colnames(geneset.clustered[[nm]])[-3] <- paste0(colnames(geneset.clustered[[nm]])[-3], " (", nm, ")")
    }
  }

  dat.xpt <- Reduce(function(...) merge(..., all=TRUE, by="gene set"), geneset.clustered)

  if (!is.null(filename)) {
    write.table(dat.xpt, filename,
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  } else {
    return(dat.xpt)
  }
}
