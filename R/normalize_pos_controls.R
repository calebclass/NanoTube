#' Positive control gene normalization
#'
#' Scale genes by the geometric mean of positive control genes.
#'
#' @param dat NanoString data, including expression matrix and gene dictionary.
#' @param logfile Optional name of logfile to print warnings or errors.

normalize_pos_controls <- function(dat, logfile="") {
  #This function scales samples by the sum of positive controls.

  # Check for positive controls
  if (sum(dat$dict$CodeClass == "Positive") == 0) {
    stop("No positive control genes found (should have CodeClass of 'Positive'). Stopping...\n",
         file=logfile)
  }

  exprs.dat <- dat$exprs

  laneGM <- apply(exprs.dat[dat$dict$CodeClass == "Positive",], 2, gm_mean)
  scale.factor <- laneGM / mean(laneGM)
  if (any(scale.factor < 0.3 | scale.factor > 3)) {
    warning(cat("Identified positive scale factor outside recommended range (0.3-3). \n
                Check samples prior to conducting analysis.\n", file=logfile, append=TRUE))
  }

  exprs.dat <- sweep(exprs.dat, 2, scale.factor, '/')
  dat$exprs <- exprs.dat
  dat$pc.scalefactors <- scale.factor
  return(dat)
}
