#' Build leading edge report from clustered mHG results
#'
#' This provides the genewise differential expression statistics and gene set membership for a
#' clustered mHG analysis.
#'
#' @param grouped.mHG The result from groupByDistance()
#' @param leadingEdge The same leading edge analysis used in groupByDistance()
#' @param limmaResults Result from runLimmaAnalysis()
#' @param filename Optional name of output text file. If NULL (default), will return as data frame.

prep.stackedReport.from.grouped.mHG <- function(grouped.mHG, leadingEdge, limmaResults,
                                                filename=NULL) {
  
  diffExpr.mat <- makeDiffExprFile(limmaResults)
  Cluster <- rep(1, times=nrow(diffExpr.mat))
  names(Cluster) <- rownames(diffExpr.mat)
  diffExpr.mat <- cbind(Cluster, diffExpr.mat)
  xpt.mat.full <- merge(diffExpr.mat, leadingEdge[,as.character(rownames(grouped.mHG))], by="row.names", all = FALSE)
  colnames(xpt.mat.full)[1] <- "Symbol"
  
  if (ncol(leadingEdge) == 1) {
    xpt.mat <- xpt.mat.full[xpt.mat.full$y == 1,]
    colnames(xpt.mat)[ncol(xpt.mat)] <- rownames(grouped.mHG)
  } else {
    #xpt.mat <- matrix(c(rep(NA, times=ncol(diffExpr.mat)+1), grouped.mHG$Cluster),
    #                  nrow=1, ncol=ncol(xpt.mat.full))
    xpt.mat <- matrix(nrow=0, ncol=ncol(xpt.mat.full))
    colnames(xpt.mat) <- colnames(xpt.mat.full)
    #colnames(xpt.mat) <- 1:ncol(xpt.mat)
    
    grouped.mHG$`gene set` <- as.character(grouped.mHG$`gene set`)
    
    for (i in 1:max(grouped.mHG$Cluster)) {
      grouped.cluster <- grouped.mHG[grouped.mHG$Cluster == i,]
      xpt.mat.full[,2] <- i
      #xpt.mat <- rbind(xpt.mat,
      #                 rbind(c(paste0("Cluster ", i), rep("", times = ncol(xpt.mat.full)-1)),
      #                 c("Symbol", colnames(diffExpr.mat), rownames(grouped.cluster),
      #                   rep("", times=sum(grouped.mHG$Cluster != i)))))
      
      leadingCluster <- xpt.mat.full[,grouped.cluster$`gene set`]
      
      if(is.null(dim(leadingCluster))) {
        xpt.mat <- rbind(xpt.mat,
                         xpt.mat.full[leadingCluster == 1,])
      } else {
        xpt.mat <- rbind(xpt.mat,
                         xpt.mat.full[rowSums(leadingCluster) >= 1,])
      }
    }
  }
  
  #xpt.mat <- xpt.mat[,colSums(xpt.mat != "") > 0]
  
  if (!is.null(filename)) {
    dir.create(filename)
    write.table(grouped.mHG, file=paste0(filename, "/summaryReport.txt"),
                col.names=TRUE, row.names=FALSE, sep = "\t", quote=FALSE)
    write.table(xpt.mat, file=paste0(filename, "/fullResults.txt"),
                col.names=TRUE, row.names=FALSE, sep = "\t", quote=FALSE)
  } else {
    return(xpt.mat)
  }
}
