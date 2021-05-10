#' Draw volcano plot of differential expression results
#'
#' Draw a volcano plot for results of a differential expression analysis by
#' limma.
#' 
#' @export
#'
#' @param limmaResults Result from runLimmaAnalysis.
#' @param plotContrast Contrast to select for volcano plot. Should be one of the
#' columns in the limma coefficients matrix (for example, a sample group that
#' was compared against the base group, or one of the contrasts in the design
#' matrix). If NULL (default), will plot the first non-Intercept column from
#' the limma coefficients matrix.
#' @param y.var The variable to plot for the y axis, either "p.value" or
#' "q.value" (the false discovery adjusted p-value)
#' 
#' @return A volcano plot using ggplot2
#' 
#' @examples
#' data(ExampleResults) # Results from runLimmaAnalysis
#' 
#' deVolcano(ExampleResults, plotContrast = "Autoimmune.retinopathy")

deVolcano <- function(limmaResults, 
                      plotContrast = NULL, y.var = c("p.value", "q.value")) {
  
    # Bind local variables
    log2FC <- log10p <- NULL
    
    # Identify contrast if not provided
    if (is.null(plotContrast)) {
        plotContrast <-
            colnames(limmaResults$coefficients)[
                which(!(colnames(limmaResults) %in% 
                        c("Intercept", "(Intercept)")))[1]]
      
        cat("\n'plotContrast' not provided, setting it to", plotContrast, "\n")
    }
    
    # Set up data frame for plot
    df <- data.frame(log2FC = limmaResults$coefficients[,plotContrast],
                     log10p = -log10(limmaResults[[y.var[1]]][,plotContrast]))
    
    plt <- ggplot(df, aes(x = log2FC, y = log10p)) +
        geom_point() +
        xlab("log2(Fold Change)") +
        ylab(paste0("-log10(", substr(y.var[1], 1, 1), ")")) +
        theme_bw()
    
    return(plt)
  
}
