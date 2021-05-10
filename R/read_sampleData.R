#' Read in a sample data table.
#'
#' Read in a .txt or .csv file containing sample names, group identifiers,
#' replicate identifiers, and any other sample data. Sample names must be in the
#' first column and must correspond with sample names in the count data file(s).
#' 
#' @export
#'
#' @param dat expression data, read in by read_merge_rcc or read.delim
#' @param file.name the path/name of the .txt or .csv file
#' @param idCol the column name of the sample identifiers in the sample table,
#' which should correspond to the column names in the count table 
#' (default NULL: will assume the first column contains the sample identifiers).
#' @param groupCol the column name of the group identifiers.
#' @param replicateCol the column name of the replicate identifiers (default
#' NULL). Multiple replicates of the same sample will have the same value in 
#' this column.
#' 
#' @return The list with the expression data, now combined with the sample
#' information
#' 
#' @examples 
#' example_data <- system.file("extdata", "GSE117751_RAW", package = "NanoTube")
#' sample_info <- system.file("extdata", "GSE117751_sample_data.csv", 
#'                            package = "NanoTube")
#' 
#' dat <- read_merge_rcc(list.files(example_data, full.names = TRUE))
#'
#' # Merge expression data with sample info
#' dat <- read_sampleData(dat, file.name = sample_info,
#'                        groupCol = "Sample_Diagnosis")

read_sampleData <- function(dat, file.name, 
                            idCol = NULL, groupCol, replicateCol = NULL) {
    # Read in the file
    file.extension <- substr(file.name, (nchar(file.name)-3), nchar(file.name))
    if (file.extension %in% c(".txt", ".TXT")) {
        sampDat <- read.delim(file.name, stringsAsFactors = FALSE)
    } else if (file.extension %in% c(".csv", ".CSV")) {
        sampDat <- read.csv(file.name, stringsAsFactors = FALSE)
    } else {
        stop("\nSample data table must be in .csv (comma-delimited) or .txt 
             (tab-delimited) format.\n")
    }
    
    # Set the sample ID's
    if (is.null(idCol)) {
        warning("\nidCol not provided. Assuming the first column of '",
                gsub(".*\\/", "", file.name), "' contains sample ID's.\n")
        
        rownames(sampDat) <- sampDat[,1]
    }
    else {
        if (idCol %in% colnames(sampDat)) {
            rownames(sampDat) <- sampDat[,idCol]
        } else {
            stop("\n'", idCol, "' is not one of the columns in '",
                 gsub(".*\\/", "", file.name), "'.\n")
        }
    }
    
    # Check if expression data and meta data have the same number of samples.
    if (ncol(dat$exprs) == nrow(sampDat)) {
        if (all(rownames(sampDat) %in% colnames(dat$exprs))) {
            # If both files have same sample names, ensure that samples 
            # are in the same order.
            dat$samples <- sampDat[colnames(dat$exprs),]
        } else {
            # If files have different sample names (or no sample names in 
            # meta file), assume that they are in the same order (and warn).
            dat$samples <- sampDat
            rownames(dat$samples) <- colnames(dat$exprs)
            
            warning("\nSample names in the two files don't match. NanoTube is 
                assuming that samples are in the same order. Please confirm 
                with your data.\n")
        }
    } else {
        stop("\nConflict: There is expression data from ", ncol(dat$exprs), 
             " samples, but sample data from ", nrow(sampDat), " samples.\n")
    }
    
    # Check that group column and replicates column are present in the table.
    if (groupCol %in% colnames(dat$samples)) {
        dat$groups <- dat$samples[,groupCol]
    } else {
        stop("\n'", groupCol, "' is not one of the column names in the
             sample table.\n")
    }
    
    if (!is.null(replicateCol)) {
        if (replicateCol %in% colnames(dat$samples)) {
            dat$replicates <- dat$samples[,replicateCol]
        } else {
            stop("\n'", replicateCol, "' is not one of the column names in the
                 sample table.\n")
        }
    }
    
    return(dat)

}