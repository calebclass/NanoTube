#' Read .gmt file.
#'
#' Read in a .gmt gene set file. From the qusage library (Bolen C).
#' 
#' @param file The filename
#'
#' @return A list object, containing a character vector of genes for
#' each gene set.
#' 
#' @examples 
#' # Write a gmt file containing 2 gene sets.
#' gene_sets <- matrix(c("Set1", "Description", "GeneA", "GeneB",
#'                       "Set2", "Description", "GeneC", "GeneD"),
#'                     nrow = 2, ncol = 4, byrow = TRUE)
#' write.table(gene_sets, file = "foo.gmt", 
#'             row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
#'             
#' # Read the gmt file
#' gene_set_list <- read_gmt("foo.gmt")
#' 
#' gene_set_list

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
