#' Postprocessing for GSEA analyses for Excel
#'
#' Clusters GSEA results by leading edge genes, and writes reports showing
#' gene expression profiles of these genes (to Excel).
#'
#' @export
#'
#' @param genesetResults Results from pathway analysis using limmaToFGSEA.
#' @param leadingEdge Results from fgseaToLEdge
#' @param limmaResults Results from runLimmaAnalysis
#' @param join.threshold The threshold distance to join gene sets. Gene sets 
#' with a distance below this value will be joined to a single "cluster."
#' @param ngroups The desired number of gene set groups. Either 
#' 'join.threshold' or 'ngroups' must be specified, 'ngroups' takes priority 
#' if both are specified.
#' @param dist.method Method for distance calculation (see options for dist()).
#' We recommend the 'binary' (also known as Jaccard) distance.
#' @param filename File name for the output Excel file.
#' 
#' @return An Excel file where the first sheet summarizes the gene set analysis 
#' results. Subsequent sheets are reports showing differential expression 
#' statistics of leading edge genes.
#' 
#' @examples
#' data("ExamplePathways")
#' data("ExampleResults") # Results from runLimmaAnalysis
#' 
#' fgseaResults <- limmaToFGSEA(ExampleResults, gene.sets = ExamplePathways)
#' 
#' leadingEdge <- fgseaToLEdge(fgseaResults, cutoff.type = "padj", cutoff = 0.1)
#' 
#' \donttest{
#' fgseaPostprocessingXLSX(fgseaResults, leadingEdge, 
#'                     limmaResults = ExampleResults,
#'                     join.threshold = 0.5,
#'                     filename = "Results.xlsx")
#' }


fgseaPostprocessingXLSX <- function(genesetResults, leadingEdge, limmaResults,
                                    join.threshold = 0.5, ngroups = NULL,
                                    dist.method = "binary", filename) {

    tab <- makeFGSEAmasterTable(genesetResults, leadingEdge,
                                  join.threshold)
    xlsx::write.xlsx(tab, filename, sheetName="Summary",
               row.names = FALSE, col.names = TRUE, append = FALSE)
    
    for (i in names(genesetResults)) {
        if(!is.null(leadingEdge[[i]])) {
            grouped <-   groupFGSEA(genesetResults[[i]],
                                         leadingEdge[[i]],
                                         join.threshold,
                                         ngroups, dist.method,
                                         returns = "signif")
            
            tab <- groupedGSEAtoStackedReport(grouped, leadingEdge[[i]], 
                                                       limmaResults)
            xlsx::write.xlsx(tab, filename, sheetName=i,
                             row.names = FALSE, col.names = TRUE, append = TRUE)
        }
    }
  
}
