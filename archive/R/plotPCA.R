#' Plot PCA
#'
#' Conduct principal components analysis and plot the results.
#'
#' @param ns Processed NanoString data
#' @param pc1 Principal component to plot on x-axis (default 1)
#' @param pc2 Principal component to plot on y-axis (default 2)
#' @param interactive.plot Plot using plotly? Default FALSE
#' @return A list containing:
#' \item{tab}{The PCA table}
#' \item{plt}{The PCA plot}

plotPCA <- function(ns, pc1 = 1, pc2 = 2, interactive.plot = FALSE) {
  pca.dat <- log2(ns$exprs[ns$dict$CodeClass == "Endogenous" &
                                 rowSums(ns$exprs == 0) == 0,]+0.5)
  
  pca <- prcomp(t(pca.dat),
                center = TRUE, scale = TRUE)
  pca.ly <- as.data.frame(pca$x[,1:2])
  pca.ly$sample <- row.names(pca$x)
  pca.ly$group <- ns$deRes$sampleData$group
  
  perc.var <- pca$sdev^2 / sum(pca$sdev^2) * 100
  
  layout.x <- layout.y <- list(
    showline = TRUE,
    showticklabels = TRUE,
    showgrid = FALSE,
    zeroline = FALSE,
    linecolor = toRGB("black"),
    linewidth = 1
  )
  
  layout.x$title <- paste0("PC1 (", signif(perc.var[1], 2), "% var)")
  layout.y$title <- paste0("PC2 (", signif(perc.var[2], 2), "% var)")
  
  plt <- plot_ly(data = pca.ly, x = ~PC1, y = ~PC2,
                 text = ~paste0("Sample: ", sample), color = ~group,
                 type = "scatter", mode = "markers") %>%
    layout(xaxis = layout.x, yaxis = layout.y)
  
  return(plt)
}