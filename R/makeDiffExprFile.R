#' Make differential expression results file.
#'
#' Make a data frame or text file containing coefficients, p-, and q-values from Limma
#' differential expression analysis. If returns == "all", will also center the log-expression
#' data on median of base.group expression and include the expression data in the output.
#'
#' @param limmaResults Result from runLimmaAnalysis
#' @param filename The desired name for the output tab-delimited text file. If NULL (default)
#' the resulting table will be returned as an R data frame.
#' @param returns If "all" (default), will center the log-expression
#' data on median of base.group expression and include the expression data in the output. If
#' "stats", will only include the differential expression statistics.

makeDiffExprFile <- function(limmaResults, filename = NULL,
                             returns = c("all", "stats"),
                             skipFirst = TRUE) {
  dat.scaled <- Biobase::exprs(limmaResults$eset) -
    apply(Biobase::exprs(limmaResults$eset)[,limmaResults$eset$group == levels(limmaResults$eset$group)[1]], 1, median)

  dat.scaled <- dat.scaled[,order(limmaResults$eset$group)]

  limma.dat <- list()

  inc <- colnames(limmaResults$t)
  if (skipFirst) inc <- inc[-1]

  for (i in inc) {
    limma.dat[[i]] <- cbind(limmaResults$coefficients[,i], limmaResults$t[,i], limmaResults$p.value[,i], limmaResults$q.value[,i])
    colnames(limma.dat[[i]]) <- paste0(c("Log2FC (", "t (", "p-val (", "q-val ("), i, ")")
  }
  limma.tab <- do.call(cbind, limma.dat)

  if (returns == "stats") {

    limma.tab <- signif(limma.tab, digits = 2)

    if (!is.null(filename)) {
      dat.xpt <- cbind(rownames(limma.tab), limma.tab)
      colnames(dat.xpt)[1] <- "Symbol"
      write.table(dat.xpt, filename,
                  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    } else {
      return(limma.tab)
    }
  } else {

    dat.xpt <- cbind(limma.tab, dat.scaled)

    if (!is.null(filename)) {
      dat.xpt <- cbind(rownames(dat.xpt), dat.xpt)
      colnames(dat.xpt)[1] <- "Symbol"
      write.table(dat.xpt, filename,
                  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    } else {
      return(dat.xpt)
    }
  }

}
