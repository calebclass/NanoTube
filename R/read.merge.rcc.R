#' Merge multiple .rcc files
#'
#' Read in multiple .rcc files named in the fileList and merge the expression data.
#'
#' @param fileList a character vector of .rcc file names
#' @param includeQC include merged QC data (from the "Lane Attributes" part of file) in
#' the output? Default FALSE
#' @param logfile a filename for the logfile (optional). If blank, will print warnings to screen.

read.merge.rcc <- function(fileList, includeQC = FALSE, logfile = "") {
  dat <- lapply(fileList, read.rcc)

  cat("Checking codeset contents......", file=logfile, append=TRUE)
  # Check that all files have same gene names in same order
  if (length(unique(lapply(dat, function(df) df$exprs$Name))) > 1) {
    stop(cat("RCC files do not have same gene names. Stopping...", file=logfile, append=TRUE))
  }
  cat("Checked gene name consistency in .RCC files.", file=logfile, append=TRUE)

  merged.dat <- list(exprs = sapply(dat, function(X) X[[1]][,4]),
                     dict = dat[[1]][[1]][,1:3])
  #                   samples = cbind(dat[[1]][[2]][,1], sapply(dat, function(X) X[[2]][,2])))
  #                     qc = cbind(dat[[1]][[3]][,1], sapply(dat, function(X) X[[3]][,2])))
  rownames(merged.dat[[1]]) <- rownames(merged.dat[[2]]) <- merged.dat$dict[,3]
  colnames(merged.dat[[1]]) <- gsub(".*\\/", "", fileList)
  #colnames(merged.dat[[3]]) <- c(colnames(dat[[3]][[1]])[1], gsub(".*\\/", "", fileList))

  if (any(merged.dat$dict$CodeClass == "Housekeeping")) {
    cat(paste(sum(merged.dat$dict$CodeClass == "Housekeeping")), "genes labeled as 'housekeeping'.",
        file=logfile, append=TRUE)
  } else {
    warning(cat("No genes labeled as 'housekeeping' in data set.", file=logfile, append=TRUE))
  }

  if (includeQC == TRUE) {
    merged.dat$qc <-  cbind(dat[[1]][[3]][,1], sapply(dat, function(X) X[[3]][,2]))
    rownames(merged.dat$qc) <- merged.dat$qc[,1]
    merged.dat$qc <- merged.dat$qc[,-1]
    colnames(merged.dat$qc) <- colnames(merged.dat$exprs)
  }
  #  colnames(merged.dat[[4]]) <- c(colnames(dat[[4]][[1]])[1], gsub(".*\\/", "", fileList))
  return(merged.dat)
}
