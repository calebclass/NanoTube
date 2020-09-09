mHG.postprocessing.xlsx <- function(genesetResults, leadingEdge, limmaResults,
                                    join.threshold = 0.5, filename) {

  tab <- make.mHG.master.table(genesetResults, leadingEdge,
                                join.threshold)
  xlsx::write.xlsx(tab, filename, sheetName="Summary",
             row.names = FALSE, col.names = TRUE, append = FALSE)
  
  for (i in names(genesetResults)) {
    for (j in names(genesetResults[[i]])) {
      if(!is.null(leadingEdge[[i]][[j]])){
        grouped <-   groupByDistance(genesetResults[[i]][[j]],
                                     leadingEdge[[i]][[j]],
                                     join.threshold)
        
        tab <- prep.stackedReport.from.grouped.mHG(grouped, leadingEdge[[i]][[j]], 
                                                   limmaResults)
        xlsx::write.xlsx(tab, filename, sheetName=paste0(i, "_", j),
                   row.names = FALSE, col.names = TRUE, append = TRUE)
      }
    }
  }

}
