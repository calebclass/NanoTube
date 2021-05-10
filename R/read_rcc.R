#' Read .rcc file
#'
#' This function reads in a single .rcc file and splits into expression, 
#' sample data, and qc components.
#' 
#' @export
#'
#' @param file file name
#' 
#' @return list containing expression data, sample attributes, and basic qc from
#' the .rcc file.
#' 
#' @examples
#' example_data <- system.file("extdata", "GSE117751_RAW", package = "NanoTube")
#' 
#' # First file only
#' single_file <- list.files(example_data, full.names = TRUE)[1]
#' single_dat <- read_rcc(single_file)

read_rcc <- function(file) {

    dat.csv <- read.csv(file,
                        row.names = NULL, col.names = seq_len(4), 
                        header = FALSE, stringsAsFactors = FALSE)
    dat <- list(exprs = dat.csv[(which(dat.csv[,1] == "<Code_Summary>")+1):
                                  (which(dat.csv[,1] == "</Code_Summary>")-1),],
                sample = dat.csv[
                    (which(dat.csv[,1] == "<Sample_Attributes>")+1):
                        (which(dat.csv[,1] == "</Sample_Attributes>")-1), 
                    c(1,2)],
                qc = dat.csv[(which(dat.csv[,1] == "<Lane_Attributes>")+1):
                               (which(dat.csv[,1] == "</Lane_Attributes>")-1), 
                             c(1,2)])
    for (i in seq_len(3)) {
        colnames(dat[[i]]) <- dat[[i]][1,]
        dat[[i]] <- dat[[i]][-1,]
        rownames(dat[[i]]) <- NULL
    }
    dat[[1]]$Count <- as.numeric(dat[[1]]$Count)
    return(dat)
}
