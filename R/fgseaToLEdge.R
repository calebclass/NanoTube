#' Generate leading edge matrix from fgsea results.
#' 
#' Extract leading edge genes from gene sets identified in fgsea analysis.
#' Gene sets may be filtered by significance or NES.
#' 
#' @export
#' 
#' @param fgsea.res Result from limmaToFGSEA
#' @param cutoff.type Filter gene sets by adjusted p-value ('padj'), nominal
#' p-value ('pval'), normalized enrichment score ('NES'), or include all gene
#' sets ('none')
#' @param cutoff Numeric cutoff for filtering (not used if 
#' cutoff.type == "none")
#' @param nes.abs.cutoff If cutoff.type == "NES", should we use extreme positive
#' and negative values (TRUE), or only filter in the positive or negative 
#' direction (FALSE). If TRUE, will select gene sets with abs(NES) > cutoff. 
#' If FALSE, will select gene sets with NES > cutoff (if cutoff >= 0) or
#' NES < cutoff (if cutoff < 0)
#' @return a list containing the leading edge matrix for each comparison
#' 
#' @examples 
#' data("ExamplePathways")
#' data("ExampleResults") # Results from runLimmaAnalysis
#' 
#' fgseaResults <- limmaToFGSEA(ExampleResults, gene.sets = ExamplePathways)
#' 
#' # Generate the leading edge for pathways with padj < 0.25
#' leadingEdge <- fgseaToLEdge(fgseaResults, 
#'                             cutoff.type = "padj", cutoff = 0.25)
#' 
#' # Generate the leading edge for pathways with abs(NES) > 2
#' leadingEdge <- fgseaToLEdge(fgseaResults, cutoff.type = "NES",
#'                             cutoff = 2, nes.abs.cutoff = TRUE)

fgseaToLEdge <- function(fgsea.res, 
                         cutoff.type = c("padj", "pval", "NES", "none"),
                         cutoff = 0.05, nes.abs.cutoff = TRUE) {
  
    # Check cutoff
    cutoff.type <- cutoff.type[1]
    if (cutoff.type %in% c("padj", "pval")) {
        if (cutoff <= 0 | cutoff > 1) {
            stop("cutoff must be between 0 and 1 for 'pval' or 'padj'.")
        }
    } else if (!(cutoff.type %in% c("padj", "pval", "NES", "none"))) {
        stop("cutoff.type must be 'padj', 'pval', 'NES', or 'none'.")
    }
    
    ledge.res <- list()
    
    for (i in names(fgsea.res)) {
        fgsea.res[[i]] <- as.data.frame(fgsea.res[[i]])
        
        if (cutoff.type %in% c("padj", "pval")) {
            fgsea.sub <- 
              fgsea.res[[i]][as.vector(fgsea.res[[i]][,cutoff.type] < cutoff & 
                                         !is.na(fgsea.res[[i]][,cutoff.type])),]
        } else if (cutoff.type == "NES") {
            if (nes.abs.cutoff) {
                fgsea.sub <- fgsea.res[[i]][abs(fgsea.res[[i]][,"NES"]) > 
                                              cutoff & 
                                            !is.na(fgsea.res[[i]][,"NES"]),]
            } else if (cutoff >= 0) {
                fgsea.sub <- fgsea.res[[i]][fgsea.res[[i]][,"NES"] > cutoff & 
                                            !is.na(fgsea.res[[i]][,"NES"]),]
            } else {
                fgsea.sub <- fgsea.res[[i]][fgsea.res[[i]][,"NES"] < cutoff & 
                                            !is.na(fgsea.res[[i]][,"NES"]),]
            }
        } else {
            fgsea.sub <- fgsea.res[[i]]
        }
        
        if (nrow(fgsea.sub) == 0) {
            # Returning a NULL value in place of a leading edge matrix
            ledge.res[[i]] <- ""
            ledge.res[i] <- list(NULL)
            warning("No gene sets meeting cutoff. Will return empty leading 
                    edge matrix")
        } else {
            ledge.dat <- fgsea.sub$leadingEdge
            ledge.genes <- unique(unlist(ledge.dat))
            
            ledge.res[[i]] <- matrix(0, nrow = length(ledge.genes), 
                                     ncol = nrow(fgsea.sub),
                                     dimnames = list(ledge.genes, 
                                                     fgsea.sub$pathway))
            
            for (j in seq_len(nrow(fgsea.sub))) 
              ledge.res[[i]][ledge.dat[[j]],j] <- 1
        }
    }
    
    return(ledge.res)
}
