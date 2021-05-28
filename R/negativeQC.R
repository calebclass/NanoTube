#' Calculate negative control statistics
#'
#' Provide a table the negative control statistics, and plot the counts of 
#' negative control genes in each sample.
#' 
#' @export
#'
#' @param ns NanoString data, processed by `processNanostringData` with
#' output.format set to 'list' and 'nSolver' normalization.
#' @param interactive.plot Generate an interactive plot using plotly? Only 
#' recommended for fewer than 20 samples (default FALSE)
#' @return A list object containing:
#' \item{tab}{The table of negative control statistics, including the mean
#' & standard deviation of negative control genes, calculated background
#' threshold, and number of endogenous genes below that threshold}
#' \item{plt}{An object containing the negative control plots.}
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
#'                              bgType = "threshold", 
#'                              bgThreshold = 2, bgProportion = 0.5,
#'                              output.format = "list")
#' 
#' negQC <- negativeQC(dat, interactive.plot = FALSE) 
#' 
#' # View negative QC table & plot
#' head(negQC$tab)
#' negQC$plt

negativeQC <- function(ns, interactive.plot = FALSE) {
  
    if (ns$normalization != "nSolver") {
        stop("Must run processNanostringData with normalization = 'nSolver' 
             prior to using this function")
    }
    
    if (!is(ns, "list")) {
        stop("Must run processNanostringData with output.format = 'list' prior
             to using this function")
    }
    
    # Bind local variables
    Count <- Sample <- Gene <- NULL
    
    # Strip plot for negative control genes
    dat.neg <- as.data.frame(ns$exprs.raw[ns$dict.raw$CodeClass == "Negative",])
    dat.neg$Gene <- ns$dict.raw$Name[ns$dict.raw$CodeClass == "Negative"]
    
    dat.neg.df <- reshape::melt(dat.neg, "Gene")
    colnames(dat.neg.df)[2:3] <- c("Sample", "Count")
    
    
    # Negative Table
    if (ncol(ns$bg.stats) == 3) {
        neg.tab <- round(ns$bg.stats, 2)
    } else {
        neg.tab <- round(ns$bg.stats[,seq_len(4)], 2)
        neg.tab$fail <- paste0(ns$bg.stats$num.less.bg, " (",
                               round(ns$bg.stats$frc.less.bg*100, 1), "%)")
        colnames(neg.tab) <- c("Mean (Neg)", "Max (Neg)", "sd (Neg)", 
                               "Background cutoff", 
                               "Genes below BG (%)")
    }
    
    neg.plot <- ggplot(data = dat.neg.df, 
                       aes(x=Count, y=Sample, 
                           text=paste0("Sample: ", Sample, "\nGene: ", Gene, 
                                       "\nCount: ", Count))) +
      geom_jitter(height = 0.2, width = 0, colour = "black", 
                  fill = "grey70", pch=21) +
      theme_classic() + ylab("") 
    
    if (interactive.plot) {
        neg.plot <- plotly::ggplotly(neg.plot, tooltip = c("text"), 
                                     width = 550, height = 400)
    }
    
    return(list(tab = neg.tab,
                plt = neg.plot))
}