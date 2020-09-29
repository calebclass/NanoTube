#' Negative control thresholding via t-test
#'
#' Conduct a one-sided t-test for each endogenous gene, comparing whether its average
#' expression across samples is significantly higher than the average expression of
#' negative control genes across samples. Genes failing this test are then removed.
#'
#' @param dat NanoString data, including expression matrix and gene dictionary.
#' @param pval p-value (from one-sided t-test) threshold to declare gene expression above
#' background expression level. Genes with p-values above this level are removed from
#' further analysis.
#' @return The updated NanoString data list object. Endogenous genes that do not pass
#' background thresholding are removed from exprs element. Additional background
#' statistics elements are also appended to this list.

remove.background.pval <- function(dat, pval){
  # This function tests whether genes are significantly higher than negative control genes
  # across samples.

  negative.mean <- apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, mean)
  negative.sd <- apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, sd)

  dat$bg.stats <- data.frame(Mean.Neg = negative.mean,
                             Max.Neg = apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, max),
                             sd.Neg = negative.sd)

  dat$gene.stats <- data.frame(row.names = dat$dict$Name,
                               t.stat = rep(NA, times=length(dat$dict$Name)),
                               p.val = rep(NA, times=length(dat$dict$Name)),
                               pass = rep(NA, times=length(dat$dict$Name)))

  rowsKeep <- ifelse(dat$dict$CodeClass == "Endogenous", yes = FALSE, no = TRUE)
  background <- c(dat$exprs[dat$dict$CodeClass == "Negative",])
  pvalvec <- c()

  for (i in which(dat$dict$CodeClass == "Endogenous")){
    ttest <- t.test(x = dat$exprs[i,], y = background, var.equal = FALSE, alternative = "greater")
    if (ttest$p.value < pval) {
      rowsKeep[i] <- TRUE
    }
    pvalvec <- c(pvalvec, ttest$p.value)
    dat$gene.stats$t.stat[i] <- as.numeric(ttest$statistic)
    dat$gene.stats$p.val[i] <- ttest$p.value
  }

  dat$gene.stats$pass <- rowsKeep
  dat$exprs <- dat$exprs[rowsKeep,]
  dat$dict <- dat$dict[rowsKeep,]
  return(dat)
}
