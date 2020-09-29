#' Read .gmt file.
#'
#' Read in a .gmt gene set file. From the qusage library.
#'
#' @param file The filename
#'
#' @return A list object, containing a character vector of genes for
#' each gene set.

read_gmt <- function(file) {
  if (!grepl("\\.gmt$", file)[1]) {
    stop("Pathway information must be a .gmt file")
  }
  geneSetDB = readLines(file)
  geneSetDB = strsplit(geneSetDB, "\t")
  names(geneSetDB) = sapply(geneSetDB, "[", 1)
  geneSetDB = lapply(geneSetDB, "[", -1:-2)
  geneSetDB = lapply(geneSetDB, function(x) {
    x[which(x != "")]
  })
  return(geneSetDB)
}
