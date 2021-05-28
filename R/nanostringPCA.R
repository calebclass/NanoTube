#' Plot PCA
#'
#' Conduct principal components analysis and plot the results, using either
#' ggplot2 or plotly.
#'
#' @export
#'
#' @param ns Processed NanoString data
#' @param pc1 Principal component to plot on x-axis (default 1)
#' @param pc2 Principal component to plot on y-axis (default 2)
#' @param interactive.plot Plot using plotly? Default FALSE (in which case
#' ggplot2 is used)
#' @param exclude.zeros Exclude genes that are not detected in all samples
#' (default TRUE)
#' @return A list containing:
#' \item{pca}{The PCA object}
#' \item{plt}{The PCA plot}
#' 
#' @examples 
#' example_data <- system.file("extdata", "GSE117751_RAW", package = "NanoTube")
#' sample_data <- system.file("extdata", "GSE117751_sample_data.csv", 
#'                            package = "NanoTube")
#' 
#' # Process and normalize data first
#' dat <- processNanostringData(example_data, 
#'                              sampleTab = sample_data, 
#'                              groupCol = "Sample_Diagnosis",
#'                              normalization = "nSolver", 
#'                              bgType = "t.test", bgPVal = 0.01)
#'                                
#' # Interactive PCA using plotly                             
#' nanostringPCA(dat, interactive.plot = TRUE)$plt
#' 
#' # Static plot using ggplot2, for the 3rd and 4th PC's.
#' nanostringPCA(dat, pc1 = 3, pc2 = 4, interactive.plot = FALSE)$plt

nanostringPCA <- function(ns, pc1 = 1, pc2 = 2, 
                          interactive.plot = FALSE, exclude.zeros = TRUE) {
  
    # Bind local variables
    PC1 <- PC2 <- group <- NULL
    
    if (is(ns, "list")) {
        pca.dat <- ns$exprs[ns$dict$CodeClass == "Endogenous",]
    } else {
        pca.dat <- exprs(ns)[fData(ns)$CodeClass == "Endogenous",]
    }
    
    if (exclude.zeros) pca.dat <- pca.dat[rowSums(pca.dat == 0) == 0,]
    
    # If RUV normalization was used, data are already log-transformed.
    # Otherwise, log-transform with a pseudocount of 0.5
    if (ns$normalization[1] != "RUV") pca.dat <- log2(pca.dat + 0.5)
    
    pca <- prcomp(t(pca.dat),
                  center = TRUE, scale = TRUE)
    
    pca.ly <- as.data.frame(pca$x[,c(pc1, pc2)])
    colnames(pca.ly) <- c("PC1", "PC2")
    pca.ly$sample <- row.names(pca$x)
    pca.ly$group <- ns$groups
    
    # Percent of variation explained by each PC
    perc.var <- pca$sdev^2 / sum(pca$sdev^2) * 100
    
    if (interactive.plot) {
      
        layout.x <- layout.y <- list(
          showline = TRUE,
          showticklabels = TRUE,
          showgrid = FALSE,
          zeroline = FALSE,
          linecolor = plotly::toRGB("black"),
          linewidth = 1
        )
        
        layout.x$title <- paste0("PC", pc1, " (", 
                                 signif(perc.var[pc1], 2), "% var)")
        layout.y$title <- paste0("PC", pc2, " (", 
                                 signif(perc.var[pc2], 2), "% var)")
        
        `%>%` <- plotly::`%>%`
        
        plt <- plotly::plot_ly(data = pca.ly, x = ~PC1, y = ~PC2,
                               text = ~paste0("Sample: ", sample), 
                               color = ~group,
                               type = "scatter", mode = "markers") %>%
          plotly::layout(xaxis = layout.x, yaxis = layout.y)
    } else {
      
        plt <- ggplot(data = pca.ly,  
                      aes(x = PC1, y = PC2, color = group)) +
          geom_point() + theme_bw() +
          xlab(paste0("PC", pc1, " (", signif(perc.var[pc1], 2), "% var)")) +
          ylab(paste0("PC", pc2, " (", signif(perc.var[pc2], 2), "% var)"))
      
    }
  
    
    return(list(pca = pca, 
                plt = plt))
}
