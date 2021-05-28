#' Run gene set enrichment analysis using DE results.
#' 
#' Use the fgsea library to run gene set enrichment analysis from the 
#' NanoStringDiff analysis results. Genes will be ranked by their log2 fold
#' changes.
#' 
#' @export
#' 
#' @param deResults Result from NanoStringDiff::glm.LRT.
#' @param gene.sets Gene set file name, in .rds (list), .gmt, or .tab format;
#' or a list object containing the gene sets. Gene names must be
#' in the same form as in the ranked.list.
#' @param sourceDB Source database to include, only if using a .tab-format 
#' geneset.file from CPDB.
#' @param min.set Number of genes required to conduct analysis on a given gene 
#' set (default = 1). If fewer than this number of genes from limmaResults are 
#' included in a gene set, that gene set will be skipped for this analysis.
#' @return A list containing data frames with the fgsea results.
#' 
#' @examples 
#' 
#' \donttest{
#'  
#' example_data <- system.file("extdata", "GSE117751_RAW", package = "NanoTube")
#' sample_data <- system.file("extdata", "GSE117751_sample_data.csv", 
#'                            package = "NanoTube")
#' 
#' datNoNorm <- processNanostringData(nsFiles = example_data,
#'                                    sampleTab = sample_data, 
#'                                    groupCol = "Sample_Diagnosis",
#'                                    normalization = "none")
#'
#' # Convert to NanoString Set, retaining 2 samples per group for this example
#' # (will run faster, but still pretty slow)
#' nsDiffSet <- makeNanoStringSetFromEset(datNoNorm[,c(1,2,15,16,29,30)])
#' 
#' # Run NanoStringDiff analysis
#' nsDiffSet <- NanoStringDiff::estNormalizationFactors(nsDiffSet)
#' result <- NanoStringDiff::glm.LRT(nsDiffSet, 
#'                                   design.full = as.matrix(pData(nsDiffSet)),
#'                                   contrast = c(1, -1, 0)) 
#'                                   #contrast: Autoimmune retinopathy vs. None
#' 
#' # FGSEA with example pathways, only for pathways with at least 5 genes
#' # analyzed in NanoString experiment
#' data("ExamplePathways")
#' fgseaResult <- nsdiffToFGSEA(result, gene.sets = ExamplePathways,
#'                              min.set = 5)
#' 
#' 
#' }

nsdiffToFGSEA <- function(deResults, gene.sets, sourceDB = NULL,
                         min.set = 1) {
  
    if (is(gene.sets, "list")) {
        gene.set.list <- gene.sets
    } else {
        file.type <- substr(gene.sets, start = nchar(gene.sets)-2, 
                            stop = nchar(gene.sets))
        if (file.type == "rds") {
            gene.set.list <- readRDS(gene.sets)
        } else if (file.type == "gmt") {
            gene.set.list <- qusage::read.gmt(gene.sets)
        } else if (file.type == "tab") {
            gene.set.list <- read_cpdb_tab(gene.sets, sourceDB)
        } else {
            stop("gene.sets must be in .rds (list object), .gmt, or 
                 .tab format")
        }
    }
  
    de.stats <- deResults$table$logFC
    names(de.stats) <- rownames(deResults$table)
    
    fgsea.res <- list(res = fgsea::fgseaMultilevel(pathways = gene.set.list,
                                                   stats = de.stats,
                                                   minSize = min.set))
    
    return(fgsea.res)
  
}
