#' Process NanoString nCounter gene expression data.
#'
#' This function reads in a zip file or folder containing multiple .rcc files, and then conducts positive
#' control normalization, background correction, and housekeeping normalization.
#'
#' @param fileDirs file directory (or zip file) containing the .rcc files, or multiple directories in
#' a character vector.
#' @param bgType type of background correction to use: "threshold" sets a thresold for N standard deviations
#' above the mean of negative controls. "t.test" conducts a one-sided t test for each gene against all
#' negative controls.
#' @param bgThreshold If bgType=="threshold", number of sd's above the mean to set as threshold
#' for background correction.
#' @param bgProportion If bgType=="threshold", proportion of samples that a gene must be above threshold to
#' be included in analysis.
#' @param bgPVal If bgType=="t.test", p-value threshold to use for gene to be included in analysis.
#' @param bgSubtract Should calculated background levels be subtracted from reported expressions?
#' If TRUE, will subtract mean+numSD*sd of the negative controls from the endogenous genes,
#' and then set negative values to zero (default FALSE)
#' @param housekeeping vector of genes (symbols or accession) to use for housekeeping correction. If NULL,
#' will use genes listed in .rcc files as "Housekeeping"
#' @param skip.housekeeping Skip housekeeping normalization? (default FALSE)
#' @param includeQC Should we include the QC when reading in .rcc files? This can cause errors,
#' particularly when reading in files from multiple experiments.
#' @param sampIds a vector of sample identifiers, important if there are technical replicates.
#' Currently, this function averages technical replicates.
#' @param logfile a filename for the logfile (optional). If blank, will print warnings to screen.
#'
#' @return An rds file containing the raw and normalized counts, sample and qc info (from rcc files), and dictionary

processNanostringData <- function(fileDirs,
                                  bgType = c("threshold", "t.test"),
                                  bgThreshold = 3, bgProportion = 0.5, bgPVal = 0.001, bgSubtract = FALSE,
                                  housekeeping = NULL, skip.housekeeping = FALSE,
                                  includeQC = FALSE,
                                  sampIds = NULL,
                                  logfile = "") {

  # Quick error checking:
  if (bgThreshold < 0) {
    warning(cat("Negative background threshold detected. This is the number of \n
                standard deviations above the background mean, and should be \n
                0 or positive. Setting to 0...\n", file=logfile, append=TRUE))
    bgThreshold <- 0
  }
  if (bgProportion < 0 | bgProportion > 1) {
    stop(cat("Proportion should be between 0 and 1, the proportion of samples \n
             that must have greater expression than background to keep for \n
             analysis. Stopping...\n", file=logfile, append=TRUE))
  }

  # Extract from fileDirs, if zipped
  if (substr(fileDirs[1], (nchar(fileDirs[1])-3), nchar(fileDirs[1])) %in% c(".zip", ".ZIP")){
    fileDirs <- unzip.dirs(fileDirs)
  }

  # Get filenames (combines files from multiple directories if necessary)
  fileNames <- c(sapply(fileDirs, list.files, full.names = TRUE))

  cat("---Running processNanostringData.R ---\nReading in .RCC files......",
      file=logfile, append=TRUE)
  # Read in .rcc files
  dat <- read.merge.rcc(fileNames, includeQC, logfile)

  # Average counts for technical replicates
  if (!is.null(sampIds) & any(duplicated(sampIds))) {
    dupSamps <- names(table(sampIds)[table(sampIds) > 1])
    for (i in dupSamps) {
      dat$exprs[,sampIds == i] <- rowMeans(dat$exprs[,sampIds == i])
    }
    # Remove duplicates
    dat$exprs <- dat$exprs[,!duplicated(sampIds)]
    dat$samples <- dat$samples[,!duplicated(sampIds)]
    if (includeQC) dat$qc <- dat$qc[,!duplicated(sampIds)]
  }

  # Mark specified genes as housekeeping (may already be marked)
  dat$dict$CodeClass[dat$dict$Name %in% housekeeping | dat$dict$Accession %in% housekeeping] <- "Housekeeping"

  # Normalize positive controls
  cat("Calculating positive scale factors......",
      file=logfile, append=TRUE)
  dat.norm <- normalize.pos.controls(dat, logfile)

  # Remove genes that fail background check
  cat("Checking endogenous genes against background threshold......",
      file=logfile, append=TRUE)
  #if (bgType == "threshold"){
  #  dat.norm <- remove.background(dat.norm, mode=bgType, numSD = bgThreshold, proportionReq = bgProportion, subtract)
  #} else if (bgType == "t.test") {
  #  dat.norm <- remove.background(dat.norm, mode=bgType, pval = bgPVal, subtract)
  #}
  dat.norm <- remove.background(dat.norm, mode = bgType, 
                                numSD = bgThreshold, proportionReq = bgProportion,
                                pval = bgPVal, subtract = bgSubtract)

  # Normalize housekeeping
  if (!skip.housekeeping) dat.norm <- normalize.housekeeping(dat.norm, housekeeping)

  # Save results
  dat.out <- dat.norm
  dat.out$exprs.raw <- dat$exprs
  dat.out$dict.raw <- dat$dict

  return(dat.out)

}
