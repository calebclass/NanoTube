#' Run gene set enrichment analysis using DE results.
#' 
#' Use the fgsea library to run gene set enrichment analysis from the 
#' NanoStringDiff analysis results. Genes will be ranked by their log2 fold
#' changes.
#' 
#' @param limmaResults Result from runLimmaAnalysis.
#' @param gene.sets Gene set file name, in .rds (list), .gmt, or .tab format;
#' or a list object containing the gene sets. Gene names must be
#' in the same form as in the ranked.list.
#' @param sourceDB Source database to include, only if using a .tab-format 
#' geneset.file from CPDB.
#' @param min.set Number of genes required to conduct analysis on a given gene 
#' set (default = 1). If fewer than this number of genes from limmaResults are 
#' included in a gene set, that gene set will be skipped for this analysis.
#' @return A list containing data frames with the fgsea results.

nsdiffToFGSEA <- function(deResults, gene.sets, sourceDB = NULL,
                         min.set = 1) {
  
  if (class(gene.sets) == "list") {
    gene.set.list <- gene.sets
  } else {
    file.type <- substr(gene.sets, start = nchar(gene.sets)-2, stop = nchar(gene.sets))
    if (file.type == "rds") {
      gene.set.list <- readRDS(gene.sets)
    } else if (file.type == "gmt") {
      gene.set.list <- read.gmt(gene.sets)
    } else if (file.type == "tab") {
      gene.set.list <- read_cpdb_tab(gene.sets, sourceDB)
    } else {
      stop("gene.sets must be in .rds (list object), .gmt, or .tab format")
    }
  }

  de.stats <- deResults$table$logFC
  names(de.stats) <- rownames(deResults$table)
  
  fgsea.res <- list(res = fgsea::fgseaMultilevel(pathways = gene.set.list,
                                                 stats = de.stats,
                                                 minSize = min.set))
  
  return(fgsea.res)
  
}
