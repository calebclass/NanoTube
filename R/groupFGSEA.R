#' Cluster gene set analysis results
#'
#' Groups the pathway analysis results (using limmaToFGSEA or nsdiffToFGSEA) 
#' based on the enriched gene sets' leading edges. If the calculated distance 
#' metric is lower than the given threshold (i.e. the gene sets have highly 
#' overlapping leading edge genes), these gene sets will be joined to a single 
#' gene set "cluster." Or if 'ngroups' is specified, gene sets will be clustered
#' by similarity into that number of groups.
#'
#' @export
#'
#' @param gsea.res Results from pathway analysis for a single comparison, 
#' using limmaToFSEA.
#' @param l.edge Leading edge result from fgseaToLEdge.
#' @param join.threshold The threshold distance to join gene sets. Gene sets 
#' with a distance below this value will be joined to a single "cluster."
#' @param ngroups The desired number of gene set groups. Either
#' 'join.threshold' or 'ngroups' must be specified, 'ngroups' takes priority
#' if both are specified.
#' @param dist.method Method for distance calculation (see options for dist()).
#' We recommend the 'binary' (also known as Jaccard) distance.
#' @param returns Either "signif" or "all". This argument defines whether only 
#' significantly enriched gene sets are included in the output table, or if 
#' the full results are included. Regardless of this selection, only 
#' significantly enriched gene sets are clustered.
#' @return A data frame including the FGSEA results, plus two additional 
#' columns for the clustering results:
#' @return \item{Cluster}{The cluster that the gene set was assigned to. Gene 
#' sets in the same cluster have a distance below the join.threshold.}
#' @return \item{best}{Whether the gene set is the most enriched (by p-value) 
#' in a given cluster.}
#'
#' @examples
#' data("ExamplePathways")
#' data("ExampleResults") # Results from runLimmaAnalysis
#'
#' fgseaResults <- limmaToFGSEA(ExampleResults, gene.sets = ExamplePathways,
#'                              min.set = 5, rank.by = "t")
#'
#' leadingEdge <- fgseaToLEdge(fgseaResults, cutoff.type = "padj", 
#'                             cutoff = 0.25)
#'
#' # Group the results, and only returns those satisfying the cutoff specified 
#' # in leadingEdge()
#' groupedResults <- groupFGSEA(fgseaResults$Autoimmune.retinopathy, 
#'                              leadingEdge$Autoimmune.retinopathy,
#'                              join.threshold = 0.5,
#'                              returns = "signif")

groupFGSEA <- function(gsea.res, l.edge,
                               join.threshold = NULL,
                               ngroups = NULL,
                               dist.method = "binary",
                               returns = c("signif", "all")) {
  
    # Convert fgsea results to be usable for this
    if (colnames(gsea.res)[2] != "p.val") {
        colnames(gsea.res)[c(2,3)] <- c("p.val", "p.adj")
        gsea.res <- as.data.frame(gsea.res[,seq_len(7)])
        rownames(gsea.res) <- gsea.res$pathway
    }
    
    gsea.sub <- gsea.res[rownames(gsea.res) %in% colnames(l.edge),]
    gsea.sub <- gsea.sub[order(gsea.sub$p.adj, gsea.sub$p.val, 
                               decreasing = FALSE),]
    l.edge <- l.edge[,rownames(gsea.sub)]
    
    d <- dist(t(l.edge), method = dist.method)
    groups <- cutree(hclust(d), h = join.threshold, k = ngroups)
    
    gsea.sub$Cluster <- groups
    gsea.sub$Cluster.Max <- ifelse(duplicated(groups), yes = "", no = "x")
    
    # Use "signif" if no returns option specified
    returns <- returns[1]
    
    if (returns == "signif") {
        gsea.return <- cbind(gsea.sub, t(l.edge))
        gsea.return <- gsea.return[order(gsea.return$Cluster, 
                                         decreasing=FALSE),]
    } else {
        gsea.non <- gsea.res[!(rownames(gsea.res) %in% colnames(l.edge)),]
        gsea.non$Cluster <- NA
        gsea.non$Cluster.Max <- ""
        gsea.return <- rbind(gsea.sub, gsea.non)
    }
  
    return(gsea.return)
    
}