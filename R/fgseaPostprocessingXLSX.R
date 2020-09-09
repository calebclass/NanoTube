

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
