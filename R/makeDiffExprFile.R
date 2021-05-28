#' Make differential expression results file.
#'
#' Make a data frame or text file containing coefficients, p-, and q-values 
#' from Limma differential expression analysis. If returns == "all", will also 
#' center the log-expression data on the median of base.group expression, and 
#' include the expression data in the output.
#' 
#' @export
#' 
#' @param limmaResults Result from runLimmaAnalysis
#' @param filename The desired name for the output tab-delimited text file. If 
#' NULL (default) the resulting table will be returned as an R data frame.
#' @param returns If "all" (default), will center the log-expression
#' data on median of base.group expression and include the expression data in 
#' the output. If "stats", will only include the differential expression 
#' statistics.
#' @param skip.first Logical: Skip the first factor for gene set analysis?
#' Frequently the first factor is the 'Intercept', which is generally 
#' uninteresting for GSEA (default TRUE).
#' 
#' @return A table of differential expression results
#' 
#' @examples 
#' data("ExampleResults") # Results from runLimmaAnalysis
#' 
#' # Include expression data in the results table
#' deResults <- makeDiffExprFile(ExampleResults, returns = "all")
#' 
#' # Only include statistics, and save to a .txt file
#' \donttest{
#' makeDiffExprFile(ExampleResults, file = "DE.txt",
#'                  returns = "stats")
#' }

makeDiffExprFile <- function(limmaResults, filename = NULL,
                             returns = c("all", "stats"),
                             skip.first = TRUE) {
  
    dat.scaled <- exprs(limmaResults$eset) -
      apply(exprs(limmaResults$eset)[,limmaResults$eset$group == 
                                       levels(limmaResults$eset$group)[1]], 
            1, median)
  
    dat.scaled <- dat.scaled[,order(limmaResults$eset$group)]
  
    limma.dat <- list()
  
    inc <- colnames(limmaResults$t)
    if (skip.first) inc <- inc[-1]
  
    for (i in inc) {
        limma.dat[[i]] <- cbind(limmaResults$coefficients[,i], 
                                limmaResults$t[,i], 
                                limmaResults$p.value[,i], 
                                limmaResults$q.value[,i])
        colnames(limma.dat[[i]]) <- 
          paste0(c("Log2FC (", "t (", "p-val (", "q-val ("), i, ")")
    }
    limma.tab <- do.call(cbind, limma.dat)
  
    if (returns == "stats") {
  
      limma.tab <- signif(limma.tab, digits = 2)
  
      if (!is.null(filename)) {
          dat.xpt <- cbind(rownames(limma.tab), limma.tab)
          colnames(dat.xpt)[1] <- "Symbol"
          write.table(dat.xpt, filename,
                      row.names = FALSE, col.names = TRUE, 
                      sep = "\t", quote = FALSE)
      } else {
          return(limma.tab)
      }
    } else {
  
        dat.xpt <- cbind(limma.tab, dat.scaled)
    
        if (!is.null(filename)) {
            dat.xpt <- cbind(rownames(dat.xpt), dat.xpt)
            colnames(dat.xpt)[1] <- "Symbol"
            write.table(dat.xpt, filename,
                        row.names = FALSE, col.names = TRUE, 
                        sep = "\t", quote = FALSE)
        } else {
            return(dat.xpt)
        }
    }

}
