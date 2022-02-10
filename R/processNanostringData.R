#' Process NanoString nCounter gene expression data.
#'
#' This function reads in a zip file or folder containing multiple .rcc files 
#' (or a txt/csv file containing raw count data), and then optionally conducts 
#' positive control normalization, background correction, and housekeeping 
#' normalization.
#' 
#' @export
#'
#' @param nsFiles file path (or zip file) containing the .rcc files, or multiple
#' directories in a character vector, or a single text/csv file containing the 
#' combined counts.
#' @param sampleTab .txt (tab-delimited) or .csv (comma-delimited) file 
#' containing sample data table (optional, default NULL)
#' @param idCol the column name of the sample identifiers in the sample table,
#' which should correspond to the column names in the count table 
#' (default NULL: will assume the first column contains the sample identifiers)
#' @param groupCol the column name of the group identifiers in the sample table
#' (required if sampleTab is provided).
#' @param replicateCol the column name of the technical replicate identifiers 
#' (default NULL). Multiple replicates of the same sample will have the same 
#' value in this column. Replicates are used to improve normalization 
#' performance in the "RUV" method; otherwise they are averaged.
#' @param normalization If "nSolver" (default), continues with background, 
#' positive control, and housekeeping control normalization steps to return
#' a NanoStringSet of normalized data. If "RUV", runs RUV normalization using 
#' controls, housekeeping genes and technical replicates. If "none", returns a 
#' NanoStringSet with the raw counts, suitable for running NanoStringDiff.
#' @param bgType Only if (normalization=="nSolver"): Type of background 
#' correction to use: "threshold" sets a threshold for N standard deviations 
#' above the mean of negative controls. "t.test" conducts a one-sided t test 
#' for each gene against all negative controls.
#' @param bgThreshold If bgType=="threshold", number of sd's above the mean to 
#' set as threshold for background correction.
#' @param bgProportion If bgType=="threshold", proportion of samples that a gene
#' must be above threshold to be included in analysis.
#' @param bgPVal If bgType=="t.test", p-value threshold to use for gene to be 
#' included in analysis.
#' @param bgSubtract Should calculated background levels be subtracted from 
#' reported expressions? If TRUE, will subtract mean+numSD*sd of the negative 
#' controls from the endogenous genes, and then set negative values to zero 
#' (default FALSE)
#' @param housekeeping vector of genes (symbols or accession) to use for 
#' housekeeping correction. If NULL, will use genes listed as "Housekeeping" 
#' under CodeClass.
#' @param skip.housekeeping Skip housekeeping normalization? (default FALSE)
#' @param includeQC Should we include the QC from the .rcc files? This can 
#' cause errors, particularly when reading in files from multiple experiments.
#' @param sampIds a vector of sample identifiers, important if there are 
#' technical replicates. Currently, this function averages technical replicates.
#' sampIds will be extracted from the replicateCol in the sampleTab, if 
#' provided.
#' @param output.format If "list", will return the normalized (optional) and raw
#' expression data, as well as various QC and relevant information tables. If 
#' "ExpressionSet" (default), will convert to an n*p ExpressionSet, with n rows
#' representing genes and p columns representing samples.
#' @param logfile a filename for the logfile (optional). If blank, will print 
#' warnings to screen.
#'
#' @return An list or ExpressionSet containing the raw and/or normalized 
#' counts, dictionary, and sample info if provided
#' 
#' @examples 
#' example_data <- system.file("extdata", "GSE117751_RAW", package = "NanoTube")
#' sample_data <- system.file("extdata", "GSE117751_sample_data.csv", 
#'                            package = "NanoTube")
#' 
#' # Process NanoString data from RCC files present in example_data folder.
#' # Use standard nCounter normalization, removing genes that do
#' # pass a t test against negative control genes with p < 0.05. Return the
#' # result as an "ExpressionSet".
#' 
#' dat <- processNanostringData(nsFiles = example_data,
#'                              sampleTab = sample_data, 
#'                              groupCol = "Sample_Diagnosis",
#'                              normalization = "nSolver",
#'                              bgType = "t.test", bgPVal = 0.01,
#'                              output.format = "ExpressionSet")
#' 
#' # Load NanoString data from a csv file (from NanoString's RCC Collector tool,
#' # for example). Skip normalization by setting 'normalization = "none"'.
#' 
#' csv_data <- system.file("extdata", "GSE117751_expression_matrix.csv", 
#'                         package = "NanoTube")
#' dat <- processNanostringData(nsFile = csv_data,
#'                               sampleTab = sample_data, 
#'                               idCol = "GEO_Accession", 
#'                               groupCol = "Sample_Diagnosis",
#'                               normalization = "none")
#'                               
#' # Load NanoString data from RCC files, using a threshold background level for
#' # removing low-expressed genes. Also, specify which genes to use for 
#' # housekeeping normalization. Save the result in "list" format (useful for
#' # some QC functions) instead of an "ExpressionSet".
#' 
#' dat <- processNanostringData(nsFiles = example_data,
#'                              sampleTab = sample_data, 
#'                              groupCol = "Sample_Diagnosis",
#'                              normalization = "nSolver",
#'                              bgType = "threshold", 
#'                              bgThreshold = 2, bgProportion = 0.5,
#'                              housekeeping = c("TUBB", "TBP", "POLR2A", 
#'                                               "GUSB", "SDHA"),
#'                              output.format = "list")

processNanostringData <- function(nsFiles,
                                  sampleTab = NULL, 
                                  idCol = NULL, groupCol = NULL, 
                                  replicateCol = NULL,
                                  normalization = c("nSolver", "RUV", "none"),
                                  bgType = c("threshold", "t.test"),
                                  bgThreshold = 2, bgProportion = 0.5, 
                                  bgPVal = 0.001, bgSubtract = FALSE,
                                  housekeeping = NULL, 
                                  skip.housekeeping = FALSE,
                                  includeQC = FALSE,
                                  sampIds = NULL,
                                  output.format = c("ExpressionSet", "list"),
                                  logfile = "") {

    # Quick error checking:
    if (bgThreshold < 0) {
        warning(cat("\nNegative background threshold detected. This is the 
                    number of \nstandard deviations above the background mean, 
                    and should be \n0 or positive. Setting to 0...\n", 
                    file=logfile, append=TRUE))
        bgThreshold <- 0
    }
    if (bgProportion < 0 | bgProportion > 1) {
        stop(cat("\nProportion should be between 0 and 1, the proportion of 
                 samples \nthat must have greater expression than background to
                 keep for \nanalysis. Stopping...\n", 
                 file=logfile, append=TRUE))
    }
  
    file.extension <- substr(nsFiles[1], 
                             (nchar(nsFiles[1])-3), nchar(nsFiles[1]))
    
    # Read in expression data from individual rcc files, 
    # or merged txt or csv files.
    
    if (file.extension %in% c(".txt", ".TXT", ".csv", ".CSV")) {
      
        # Read in merged count data
        cat("\nLoading count data......",
            file=logfile, append=TRUE)
        
        if (file.extension %in% c(".txt", ".TXT")) {
            tabData <- read.delim(nsFiles,
                                stringsAsFactors = FALSE)
        } else {
            tabData <- read.csv(nsFiles,
                                stringsAsFactors = FALSE)
        }
        
        dat <- list(exprs = tabData[,-seq_len(3)],
                    dict = tabData[,seq_len(3)])
        rownames(dat$exprs) <- rownames(dat$dict) <- tabData$Accession
        
        # Remove periods or spaces from dictionary colnames
        colnames(dat$dict) <- gsub("\\.| ", "", colnames(dat$dict))
      
    } else {
      
        # Extract from nsFiles, if zipped
        if (file.extension %in% c(".zip", ".ZIP")){
            nsFiles <- unzip_dirs(nsFiles)
        }
        
        # Get filenames (combines files from multiple directories if necessary)
        fileNames <- unlist(lapply(nsFiles, list.files, full.names = TRUE))
        
        # Retain only fileNames ending in "RCC"
        # ('tolower' function makes it case-insensitive)
        fileNames <- fileNames[
          tolower(substr(fileNames, 
                         start = nchar(fileNames)-2, 
                         stop = nchar(fileNames))) == "rcc"]
        
        cat("\nReading in .RCC files......", file=logfile, append=TRUE)
        
        # Read in RCC files
        dat <- read_merge_rcc(fileNames, includeQC, logfile)
    }
    
    # Read in sample data file, if provided.
    if (!is.null(sampleTab)) dat <- read_sampleData(dat, file.name = sampleTab,
                                      idCol = idCol, groupCol = groupCol, 
                                      replicateCol = replicateCol)
    
    # Mark specified genes as housekeeping (may already be marked)
    dat$dict$CodeClass[dat$dict$Name %in% housekeeping | 
                         dat$dict$Accession %in% housekeeping] <- "Housekeeping"
    
    # Average counts for technical replicates or normalize using RUV.
    # Use the replicate ID's extracted from the sampleTab, if applicable.
    if ("replicates" %in% names(dat)) sampIds <- dat$replicates
    
    # Or use the ID's provided directly.
    # Then average replicates, or normalize using "RUV".
    normalization <- normalization[1]
    bgType <- bgType[1]
    output.format <- output.format[1]
    
    if ((!is.null(sampIds) & any(duplicated(sampIds))) |
        normalization != "RUV") {
      if (is.null(sampIds)) sampIds <- seq_len(ncol(dat$exprs))
        
        dupSamps <- names(table(sampIds)[table(sampIds) > 1])
        
        if (normalization == "RUV") {
            cat("\nConducting RUV normalization......",
                file=logfile, append=TRUE)
            
            # Save a copy of the raw expression data
            dat$exprs.raw <- dat$exprs
            dat$samples.raw <- dat$samples
            dat$dict.raw <- dat$dict
            
            # Normalize using RUV: Use all genes as control genes to start 
            # (recommended by RUV authors)
            dat$exprs <- t(ruv::RUVIII(t(log2(dat$exprs+0.5)), M = sampIds,
                                       ctl = seq_len(nrow(dat$exprs)), 
                                       k = NULL))
          
        } else {
            cat("\nAveraging technical replicates.....",
                file=logfile, append=TRUE)
            for (i in dupSamps) {
                dat$exprs[,sampIds == i] <- rowMeans(dat$exprs[,sampIds == i])
            }
        }
        
        # Remove duplicates
        dat$exprs <- dat$exprs[,!duplicated(sampIds)]
        dat$samples <- dat$samples[!duplicated(sampIds),]
        if (includeQC) dat$qc <- dat$qc[!duplicated(sampIds),]
        if ("groups" %in% names(dat))
            dat$groups <- dat$groups[!duplicated(sampIds)]
    }
    
    # Normalize using nSolver recommended method:
    if (normalization == "nSolver") {
        cat("\nCalculating positive scale factors......",
            file=logfile, append=TRUE)
        dat.norm <- normalize_pos_controls(dat, logfile)
        
        # Remove genes that fail background check
        cat("\nChecking endogenous genes against background threshold......",
            file=logfile, append=TRUE)
    
        dat.norm <- remove_background(dat.norm, mode = bgType, 
                                      numSD = bgThreshold, 
                                      proportionReq = bgProportion,
                                      pval = bgPVal, subtract = bgSubtract)
        
        if (!skip.housekeeping) {
            cat("\nConducting housekeeping normalization......",
                file=logfile, append=TRUE)
            dat.norm <- normalize_housekeeping(dat.norm, housekeeping)
        }
        
        # Save a copy of the raw data (only used if output.format == "list")
        dat.out <- dat.norm
        dat.out$exprs.raw <- dat$exprs
        dat.out$dict.raw <- dat$dict
        dat <- dat.out
    }
    
    if (output.format == "list") {
        dat$normalization <- normalization
        
        if (normalization == "none") {
            dat$exprs.raw <- dat$exprs
            dat$dict.raw <- dat$dict
        }
        
        return(dat)
      
    } else {
        dat.out <- ExpressionSet(assayData = as.matrix(dat$exprs),
                                 featureData = AnnotatedDataFrame(dat$dict))
        
        if ("samples" %in% names(dat))
            phenoData(dat.out) <- AnnotatedDataFrame(dat$samples)
        if ("groups" %in% names(dat)) phenoData(dat.out)$groups <- dat$groups
        
        phenoData(dat.out)$normalization <- normalization
        
        if (includeQC)
          phenoData(dat.out) <-
              AnnotatedDataFrame(cbind(pData(dat.out), dat$qc))
        
        return(dat.out)
    }
}
