#' Housekeeping gene normalization
#'
#' Scale endogenous genes by the geometric mean of housekeeping genes. This 
#' should be conducted after positive control normalization and background 
#' correction.  This step is conducted within processNanostringData, when 
#' normalization is set to "nCounter".
#'
#' @export
#' 
#' @param dat NanoString data, including expression matrix and gene dictionary.
#' @param genes List of housekeeping genes to use for normalization. If NULL 
#' (default), will use all genes marked as "Housekeeping" in codeset.
#' @param logfile Optional name of logfile to print messages, warnings 
#' or errors.
#' 
#' @return NanoString data, with expression matrix now normalized by 
#' housekeeping gene expression.
#' 
#' @examples
#' example_data <- system.file("extdata", "GSE117751_RAW", package = "NanoTube")
#' 
#' # Load data, positive control normalization, and background filtering
#' dat <- read_merge_rcc(list.files(example_data, full.names = TRUE))
#' dat <- normalize_pos_controls(dat)
#' dat <- remove_background(dat, mode = "t.test", pval = 0.05)
#' 
#' # Normalize by genes marked "Housekeeping" in RCC files
#' dat <- normalize_housekeeping(dat)
#' 
#' # Normalize by specified housekeeping genes (gene symbol or accession)
#' dat <- normalize_housekeeping(dat,
#'                        genes = c("TUBB", "TBP", "POLR2A", "GUSB", "SDHA"))

normalize_housekeeping <- function(dat, genes = NULL, logfile = "") {

    if (is.null(genes)) genes <- 
        dat$dict$Name[dat$dict$CodeClass == "Housekeeping"]
  
    exprs.dat <- dat$exprs
  
    if (all(genes %in% dat$dict$Name)) {
        laneGM <- apply(exprs.dat[dat$dict$Name %in% genes,], 2, gm_mean)
    } else if (all(genes %in% dat$dict$Accession)) {
        laneGM <- apply(exprs.dat[dat$dict$Accession %in% genes,], 2, gm_mean)
    } else {
        stop(cat("\nHousekeeping gene(s) not found in dataset:", 
                 genes[!(genes %in% dat$dict$Name | 
                           genes %in% dat$dict$Accession)],
                 "\n", file=logfile, append=TRUE))
    }
  
    scale.factor <- laneGM / mean(laneGM)
    exprs.dat[dat$dict$CodeClass == "Endogenous",] <-
      sweep(exprs.dat[dat$dict$CodeClass == "Endogenous",], 
            2, scale.factor, '/')
    dat$exprs <- exprs.dat
    dat$hk.scalefactors <- scale.factor
    return(dat)
}
