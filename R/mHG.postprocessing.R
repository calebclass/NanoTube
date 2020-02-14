mHG.postprocessing <- function(genesetResults, leadingEdge, limmaResults,
                               join.threshold = 0.5, reportDir) {

  dir.create(reportDir, showWarnings = FALSE)

  make.mHG.master.table(genesetResults, leadingEdge,
                        join.threshold, filename = paste0(reportDir, "/mHG_summary.txt"))

  for (i in names(genesetResults)) {
    for (j in names(genesetResults[[i]])) {
      if(!is.null(leadingEdge[[i]][[j]])){
        grouped <-   groupByDistance(genesetResults[[i]][[j]],
                                     leadingEdge[[i]][[j]],
                                     join.threshold)
        prep.stackedReport.from.grouped.mHG(grouped, leadingEdge[[i]][[j]], limmaResults,
                                            filename =  paste0(reportDir, "/", i, "_", j))
      }
    }
  }

}
