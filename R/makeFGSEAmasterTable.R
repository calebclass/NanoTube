#' Make master table of all GSEA results
#'
#' This function clusters GSEA results by leading edge similarity, and then 
#' combines to a data frame or text file.
#'
#' @param genesetResults Results from pathway analysis using limmaToFGSEA.
#' @param leadingEdge Results from fgseaToLEdge
#' @param join.threshold The threshold distance to join gene sets. Gene sets 
#' with a distance below this value will be joined to a single "cluster."
#' @param ngroups The desired number of gene set groups. Either 
#' 'join.threshold' or 'ngroups' must be specified, 'ngroups' takes priority 
#' if both are specified.
#' @param dist.method Method for distance calculation (see options for dist()).
#' We recommend the 'binary' (also known as Jaccard) distance.
#' @param filename File name for the output text file. If NULL (default), data 
#' will be returned as an R data frame.
#' 
#' @return A table of GSEA results, clustered by similarity of leading edge.

makeFGSEAmasterTable <- function(genesetResults, leadingEdge,
                                 join.threshold = 0.5,
                                 ngroups = NULL,
                                 dist.method = "binary",
                                 filename = NULL) {

    geneset.clustered <- list()
  
    for (i in names(genesetResults)) {
      
        if (is.null(leadingEdge[[i]])) {
            temp <-  genesetResults[[i]]
            temp$Cluster <- NA
            temp$Cluster.Max <- ""
            colnames(temp)[2:3] <- c("p.val", "p.adj")
            geneset.clustered[[i]] <- temp[,-8]
        } else {
            geneset.clustered[[i]] <- groupFGSEA(genesetResults[[i]], 
                                                 leadingEdge[[i]],
                                                 join.threshold, ngroups,
                                                 dist.method, returns="all")
        }
        
        colnames(geneset.clustered[[i]])[-1] <- 
          paste0(colnames(geneset.clustered[[i]])[-1], " (", i, ")")
    }
  
    dat.xpt <- Reduce(function(...) merge(..., all=TRUE, by="pathway"), 
                      geneset.clustered)
  
    if (!is.null(filename)) {
        write.table(dat.xpt, filename,
                    row.names = FALSE, col.names = TRUE, 
                    sep = "\t", quote = FALSE)
    } else {
        return(dat.xpt)
    }
}
