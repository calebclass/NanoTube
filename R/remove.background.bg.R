#' Assess background expression
#'
#' Compare endogenous gene expression data against negative control genes.
#'
#' @param dat Positive control-scaled NanoString data
#' @param numSD Number of standard deviations above mean of negative control genes to used as
#' background threshold for each sample: mean(negative_controls) + numSD * sd(negative_controls)
#' @param proportionReq Required proportion of sample expressions exceeding the sample background
#' threshold to include gene in further analysis
#' @param subtract Should calculated background levels be subtracted from reported expressions?
#' If TRUE, will subtract and then set negative values to zero (default FALSE)

remove.background.bg <- function(dat, numSD, proportionReq,
                              subtract = FALSE) {
  # This function removes genes where less than the required proportion of sample expressions exceed
  # the background threshold: mean(negative_controls) + numSD * sd(negative_controls)

  negative.mean <- apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, mean)
  negative.sd <- apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, sd)
  bg.threshold <- negative.mean + numSD * negative.sd

  dat$bg.stats <- data.frame(Mean.Neg = negative.mean,
                             Max.Neg = apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, max),
                             sd.Neg = negative.sd,
                             background = bg.threshold,
                             num.less.bg = colSums(t(t(dat$exprs[dat$dict$CodeClass == "Endogenous" ,]) < bg.threshold)),
                             frc.less.bg = colMeans(t(t(dat$exprs[dat$dict$CodeClass == "Endogenous" ,]) < bg.threshold)))

  if (subtract) {
    for (i in 1:ncol(dat$exprs)) dat$exprs[,i] <- dat$exprs[,i] - bg.threshold[i]
    dat$exprs[dat$exprs < 0] <- 0
  }

  rowsKeep <- which(dat$dict$CodeClass != "Endogenous" | rowMeans(t(t(dat$exprs) > bg.threshold)) >= proportionReq)
  dat$exprs <- dat$exprs[rowsKeep,]
  dat$dict <- dat$dict[rowsKeep,]
  return(dat)
}
