#' Build a report from gene set enrichment results.
#' 
#' After clustering FGSEA results by gene set similarity, this function builds
#' a report containing the individual gene expression profiles for genes 
#' contained in each gene set cluster.
#' 
#' @export
#' 
#' @param grouped.gsea Output from groupFGSEA()
#' @param leadingEdge Leading edge analysis results used in groupFGSEA()
#' @param de.fit Differential Expression results from Limma or NanoStringDiff
#' @param outputDir Directory for output files. If NULL (default), will return
#' the stacked report instead of writing to a file.
#' 
#' @return A stacked report containing statistics and gene expression profiles
#' for genes contained in each cluster
#' 
#' @examples
#' data("ExamplePathways")
#' data("ExampleResults") # Results from runLimmaAnalysis
#' 
#' fgseaResults <- limmaToFGSEA(ExampleResults, gene.sets = ExamplePathways,
#'                              min.set = 5, rank.by = "t")
#' leadingEdge <- fgseaToLEdge(fgseaResults, cutoff.type = "padj", cutoff = 0.1)
#' 
#' fgseaGrouped <- groupFGSEA(fgseaResults$Autoimmune.retinopathy, 
#'                             leadingEdge$Autoimmune.retinopathy,
#'                             join.threshold = 0.5,
#'                             dist.method = "binary")
#' 
#' results.AR <- groupedGSEAtoStackedReport(
#'               fgseaGrouped,
#'               leadingEdge = leadingEdge$Autoimmune.retinopathy,
#'               de.fit = ExampleResults)

groupedGSEAtoStackedReport <- function(grouped.gsea, leadingEdge, de.fit,
                                       outputDir=NULL) {
  
    diffExpr.mat <- makeDiffExprFile(de.fit, returns = "all")
    Cluster <- rep(1, times = nrow(diffExpr.mat))
    names(Cluster) <- rownames(diffExpr.mat)
    diffExpr.mat <- cbind(Cluster, diffExpr.mat)
    
    grouped.gsea <-
      grouped.gsea[rownames(grouped.gsea) %in% colnames(leadingEdge), ]
    
    xpt.mat.full <- merge(diffExpr.mat,
                          leadingEdge[, as.character(rownames(grouped.gsea))],
                          by = "row.names", all = FALSE)
    colnames(xpt.mat.full)[1] <- "Symbol"
    
    if (ncol(leadingEdge) == 1) {
        xpt.mat <- xpt.mat.full[xpt.mat.full$y == 1,]
        colnames(xpt.mat)[ncol(xpt.mat)] <- rownames(grouped.gsea)
    } else {
        xpt.mat <- matrix(nrow=0, ncol=ncol(xpt.mat.full))
        colnames(xpt.mat) <- colnames(xpt.mat.full)
    
        grouped.gsea <- cbind(as.character(rownames(grouped.gsea)), 
                              grouped.gsea)
        colnames(grouped.gsea)[1] <- 'gene set'
      
      for (i in seq_len(max(grouped.gsea$Cluster))) {
          grouped.cluster <- grouped.gsea[grouped.gsea$Cluster == i,]
          xpt.mat.full[,2] <- i
          
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
    
    if (!is.null(outputDir)) {
        dir.create(outputDir)
        write.table(grouped.gsea, file=paste0(outputDir, "/summaryReport.txt"),
                    col.names=TRUE, row.names=FALSE, sep = "\t", quote=FALSE)
        write.table(xpt.mat, file=paste0(outputDir, "/fullResults.txt"),
                    col.names=TRUE, row.names=FALSE, sep = "\t", quote=FALSE)
    } else {
        return(xpt.mat)
    }
}
