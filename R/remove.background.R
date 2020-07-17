#' Assess background expression
#'
#' Compare endogenous gene expression data against negative control genes.
#'
#' @param dat Positive control-scaled NanoString data
#' @param mode Either "threshold" (default) or "t.test". If "threshold", requires proportionReq
#' of samples to have expression numSD standard deviations among the mean of negative control
#' genes. If "t.test", 
#' @param numSD Number of standard deviations above mean of negative control genes to used as
#' background threshold for each sample: mean(negative_controls) + numSD * sd(negative_controls).
#' Required if mode == "threshold" or subtract == TRUE
#' @param proportionReq Required proportion of sample expressions exceeding the sample background
#' threshold to include gene in further analysis. Required if mode == "threshold" or 
#' subtract == TRUE
#' @param pval p-value (from one-sided t-test) threshold to declare gene expression above
#' background expression level. Genes with p-values above this level are removed from
#' further analysis. Required if mode == "t.test"
#' @param subtract Should calculated background levels be subtracted from reported expressions?
#' If TRUE, will subtract mean+numSD*sd of the negative controls from the endogenous genes,
#' and then set negative values to zero (default FALSE)

remove.background <- function(dat, 
                              mode = c("threshold", "t.test"), 
                              numSD, proportionReq, pval,
                              subtract = FALSE) {
  # This function removes genes where less than the required proportion of sample expressions exceed
  # the background threshold: mean(negative_controls) + numSD * sd(negative_controls)
  
  negative.mean <- apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, mean)
  negative.sd <- apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, sd)
  
  if (mode == "threshold" | subtract) {
    bg.threshold <- negative.mean + numSD * negative.sd
    
    dat$bg.stats <- data.frame(Mean.Neg = negative.mean,
                               Max.Neg = apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, max),
                               sd.Neg = negative.sd,
                               background = bg.threshold,
                               num.less.bg = colSums(t(t(dat$exprs[dat$dict$CodeClass == "Endogenous" ,]) < bg.threshold)),
                               frc.less.bg = colMeans(t(t(dat$exprs[dat$dict$CodeClass == "Endogenous" ,]) < bg.threshold)))
  } else {
    dat$bg.stats <- data.frame(Mean.Neg = negative.mean,
                               Max.Neg = apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, max),
                               sd.Neg = negative.sd)
  }
  
  if (mode == "t.test") {
    dat$gene.stats <- data.frame(row.names = dat$dict$Name,
                                 t.stat = rep(NA, times=length(dat$dict$Name)),
                                 p.val = rep(NA, times=length(dat$dict$Name)),
                                 pass = rep(NA, times=length(dat$dict$Name)))
    
    rowsKeep <- ifelse(dat$dict$CodeClass == "Endogenous", yes = FALSE, no = TRUE)
    background <- unlist(dat$exprs[dat$dict$CodeClass == "Negative",])
#    pvalvec <- c()
    
    for (i in which(dat$dict$CodeClass == "Endogenous")){
      ttest <- t.test(x = dat$exprs[i,], y = background, var.equal = FALSE, alternative = "greater")
      if (ttest$p.value < pval) {
        rowsKeep[i] <- TRUE
      }
#      pvalvec <- c(pvalvec, ttest$p.value)
      dat$gene.stats$t.stat[i] <- as.numeric(ttest$statistic)
      dat$gene.stats$p.val[i] <- ttest$p.value
    }
    
    dat$gene.stats$pass <- rowsKeep
  } else {
    rowsKeep <- which(dat$dict$CodeClass != "Endogenous" | 
                        rowMeans(t(t(dat$exprs) > bg.threshold)) >= proportionReq)
  }
  
  if (subtract) {
    for (i in 1:ncol(dat$exprs)) dat$exprs[,i] <- dat$exprs[,i] - bg.threshold[i]
    dat$exprs[dat$exprs < 0] <- 0
  }
  
  dat$exprs <- dat$exprs[rowsKeep,]
  dat$dict <- dat$dict[rowsKeep,]
  return(dat)
}
