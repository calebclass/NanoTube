#' Read .tab file.
#'
#' Read in a .tab file from the Consensus Pathway Database (CPDB)
#'
#' @param file The filename
#' @param sourceDB The source database to use. If NULL (default), retains
#' gene sets from all source databases
#'
#' @return A list object, containing a character vector of genes for
#' each gene set.

read_cpdb_tab <- function(file, sourceDB = NULL) {
    if (!grepl("\\.tab$", file)[1]) {
        stop("Pathway information must be a .tab file")
    }
    
    geneSetTable <- read.delim(file,
                               stringsAsFactors = FALSE)
    
    if (!is.null(sourceDB)) {
        if (!(sourceDB %in% geneSetTable$source)) {
            stop("'", sourceDB, 
                 "' is not one of the source databases in ", 
                 file)
        }
      
        geneSetTable <- geneSetTable[geneSetTable$source == sourceDB,]
    }
    
    geneSetDB <- strsplit(geneSetTable[,ncol(geneSetTable)], split = ",")
    
    if (is.null(sourceDB)) {
        names(geneSetDB) <- paste0(geneSetTable$pathway, " | ", geneSetTable$source)
    } else {
        names(geneSetDB) <- geneSetTable$pathway
    }
    
    return(geneSetDB)
}