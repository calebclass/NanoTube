#' Identify source databases from a .tab file
#'
#' Read in a .tab file from the Consensus Pathway Database (CPDB), and
#' identify the source databases present.
#'
#' @param file The filename
#'
#' @return A table of the source databases, with the number of gene sets from
#' each one.

read_cpdb_sourceDBs <- function(file) {
    if (!grepl("\\.tab$", file)[1]) {
        stop("Pathway information must be a .tab file")
    }
    
    geneSetTable <- read.delim(file,
                               stringsAsFactors = FALSE)
    
    return(table(geneSetTable$source))
}