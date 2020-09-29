#' Cluster gene set analysis results
#'
#' Groups the pathway analysis results (using mHG.genesets or mHG.genesets.multi) based on
#' the enriched gene sets' leading edges. If the calculated distance metric is lower than
#' the given threshold (i.e. the gene sets have highly overlapping leading edge genes),
#' these gene sets will be joined to a single gene set "cluster."
#'
#' @param mHG.res Results from pathway analysis for a single comparison, using mHG.genesets
#' or one member of mHG.geneset.multi. If an mHG.genesets.multi output, a single comparison
#' from the output list must be selected.
#' @param l.edge Leading edge result from mHG.to.LEdge or a single member of mHG.to.LEdge.multi.
#' If an mHG.to.LEdge.multi output, the matching comparison from mHG.res must be selected.
#' @param join.threshold The threshold distance to join gene sets. Gene sets with a distance
#' below this value will be joined to a single "cluster."
#' @param returns Either "signif" or "all". This argument defines whether only significantly
#' enriched gene sets are included in the output table, or if the full results are included.
#' Regardless of this selection, only significantly enriched gene sets (as defined in
#' mHG.to.LEdge or mHG.to.LEdge.multi) are clustered.
#' @return A data frame including the mHG results, plus two additional columns for the clustering
#' results
#' @return \item{Cluster} The cluster that the gene set was assigned to. Gene sets in the same
#' cluster have a distance below the join.threshold.
#' @return \item{best} Whether the gene set is the most enriched (by p-value) in a given cluster.

groupByDistance <- function(mHG.res, l.edge,
                            join.threshold = NULL,
                            dist.method = "binary",
                            returns = c("signif", "all")) {

  mHG.sub <- mHG.res[rownames(mHG.res) %in% colnames(l.edge),]

  if (nrow(mHG.sub) >= 2) {
    mHG.sub <- mHG.sub[order(mHG.sub$p.adj, mHG.sub$p.val, decreasing = FALSE),]
    l.edge <- l.edge[,rownames(mHG.sub)]

    d <- dist(t(l.edge), method = dist.method)
    groups <- cutree(hclust(d), h = join.threshold)

    mHG.sub$Cluster <- groups
    mHG.sub$Cluster.Max <- ifelse(duplicated(groups), yes = "", no = "x")

    mHG.sub$p.val <- signif(mHG.sub$p.val, digits=2)
    mHG.sub$p.adj <- signif(mHG.sub$p.adj, digits=2)

    mHG.return <- mHG.sub[order(mHG.sub$Cluster, decreasing=FALSE),-1]
  } else if (nrow(mHG.sub) == 1) {
    mHG.return <- mHG.sub[,-1]
    mHG.return$Cluster <- 1
    mHG.return$Cluster.Max <- "x"
  }


  #rownames(mHG.return) <- NULL

  if (returns == "all"){
    mHG.non <- mHG.res[!(rownames(mHG.res) %in% colnames(l.edge)),-1]
    mHG.non$Cluster <- NA
    mHG.non$Cluster.Max <- ""
    mHG.return <- rbind(mHG.return, mHG.non)
  }

  mHG.return <- cbind(rownames(mHG.return), mHG.return)
  colnames(mHG.return)[c(1:4,8)] <- c("gene set","genes", "index", "count", "best")
  return(mHG.return[,c(7,8,1,5,6,2:4)])

}
