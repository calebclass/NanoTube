#' Positive control gene normalization
#'
#' Scale genes by the geometric mean of positive control genes. This step is
#' conducted within processNanostringData, when normalization is set to 
#' "nCounter".
#' 
#' @export
#'
#' @param dat NanoString data, including expression matrix and gene dictionary.
#' @param logfile Optional name of logfile to print messages, warnings or 
#' errors.
#' 
#' @return NanoString data, with expression matrix now normalized by 
#' positive control gene expression.
#' 
#' @examples
#' example_data <- system.file("extdata", "GSE117751_RAW", package = "NanoTube")
#' 
#' dat <- read_merge_rcc(list.files(example_data, full.names = TRUE))
#' 
#' # Positive controls are identified in the RCC files, and used to 
#' # normalize the data
#' dat <- normalize_pos_controls(dat)

normalize_pos_controls <- function(dat, logfile="") {

    # Check for positive controls
    if (sum(dat$dict$CodeClass == "Positive") == 0) {
        stop("No positive control genes found (should have CodeClass of 
             'Positive'). Stopping...\n",
             file=logfile)
    }
  
    exprs.dat <- dat$exprs
  
    laneGM <- apply(exprs.dat[dat$dict$CodeClass == "Positive",], 2, gm_mean)
    scale.factor <- laneGM / mean(laneGM)
    
    if (any(scale.factor < 0.3 | scale.factor > 3)) {
        warning(cat("Identified positive scale factor outside recommended range 
                    (0.3-3). \nCheck samples prior to conducting analysis.\n", 
                    file=logfile, append=TRUE))
    }
  
    exprs.dat <- sweep(exprs.dat, 2, scale.factor, '/')
    dat$exprs <- exprs.dat
    dat$pc.scalefactors <- scale.factor
    return(dat)
}
