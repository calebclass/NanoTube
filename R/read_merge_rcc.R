#' Merge multiple .rcc files
#'
#' Read in multiple .rcc files named in the fileList and merge the expression 
#' data. This step is conducted within processNanostringData.
#' 
#' @export
#'
#' @param fileList a character vector of .rcc file names
#' @param includeQC include merged QC data (from the "Lane Attributes" part of 
#' file) in the output? Default FALSE
#' @param logfile a filename for the logfile (optional). If blank, will print 
#' warnings to screen.
#' 
#' @return A list object including:
#' \item{exprs}{The expression matrix}
#' \item{dict}{The gene dictionary}
#' \item{qc}{QC metrics included in the .rcc files, if includeQC == TRUE}
#' 
#' @examples 
#' example_data <- system.file("extdata", "GSE117751_RAW", package = "NanoTube")
#' 
#' dat <- read_merge_rcc(list.files(example_data, full.names = TRUE))

read_merge_rcc <- function(fileList, includeQC = FALSE, logfile = "") {
    dat <- lapply(fileList, read_rcc)
  
    cat("\nChecking codeset contents......", file=logfile, append=TRUE)
    
    # Check that all files have same gene names in same order
    # We could allow merging of discordant dictionaries in a future update.
    # Otherwise, these can be read separately and then merged.
    if (length(unique(lapply(dat, function(df) df$exprs$Name))) > 1) {
        stop(cat("\nRCC files do not have same gene names. Stopping...", 
                 file=logfile, append=TRUE))
    }
    cat("\nChecked gene name consistency in .RCC files.", 
        file=logfile, append=TRUE)
  
    vec.template <- dat[[1]][[1]][,4]
    merged.dat <- list(exprs = vapply(dat, function(X) X[[1]][,4],
                                      FUN.VALUE = vec.template),
                       dict = dat[[1]][[1]][,seq_len(3)])
  
    # We use Accessions for row names, but it's possible for those to be 
    # duplicated. To account for this, use the `make.names` function.
    rownames(merged.dat[[1]]) <- rownames(merged.dat[[2]]) <- 
      make.names(merged.dat$dict[,3], unique = TRUE)
    colnames(merged.dat[[1]]) <- gsub(".*\\/", "", fileList)
  
    if (any(merged.dat$dict$CodeClass == "Housekeeping")) {
        cat("\n", paste(sum(merged.dat$dict$CodeClass == "Housekeeping")), 
            " genes labeled as 'housekeeping'.",
            file=logfile, sep = "", append = TRUE)
    } else {
        warning(cat("\nNo genes labeled as 'housekeeping' in data set.", 
                    file=logfile, append=TRUE))
    }
  
    if (includeQC == TRUE) {
        vec.template <- dat[[1]][[3]][,2]
        merged.dat$qc <-  cbind(dat[[1]][[3]][,1], 
                                vapply(dat, function(X) X[[3]][,2],
                                       FUN.VALUE = vec.template))
        rownames(merged.dat$qc) <- merged.dat$qc[,1]
        merged.dat$qc <- t(merged.dat$qc[,-1])
        rownames(merged.dat$qc) <- colnames(merged.dat$exprs)
    }
    return(merged.dat)
}
