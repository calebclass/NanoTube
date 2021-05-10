#' Assess background expression
#'
#' Compare endogenous gene expression data against negative control genes and
#' remove data for genes that fail the comparison. This step is
#' conducted within processNanostringData, when normalization is set to 
#' "nCounter".
#' 
#' @export
#'
#' @param dat Positive control-scaled NanoString data
#' @param mode Either "threshold" (default) or "t.test". If "threshold", 
#' requires proportionReq of samples to have expression numSD standard 
#' deviations among the mean of negative control genes. If "t.test", each gene 
#' will be compared with all negative control genes in a one-sided two-sample 
#' t-test. 
#' @param numSD Number of standard deviations above mean of negative control 
#' genes to used as background threshold for each sample: 
#' mean(negative_controls) + numSD * sd(negative_controls). 
#' Required if mode == "threshold" or subtract == TRUE
#' @param proportionReq Required proportion of sample expressions exceeding the
#' sample background threshold to include gene in further analysis. Required if
#' mode == "threshold" or subtract == TRUE
#' @param pval p-value (from one-sided t-test) threshold to declare gene 
#' expression above background expression level. Genes with p-values above 
#' this level are removed from further analysis. Required if mode == "t.test"
#' @param subtract Should calculated background levels be subtracted from 
#' reported expressions? If TRUE, will subtract mean+numSD*sd of the negative 
#' controls from the endogenous genes, and then set negative values to zero 
#' (default FALSE).
#' 
#' @return NanoString data, with genes removed that fail the comparison test
#' against negative control genes. Expression levels are updated for all genes
#' if subtract == TRUE.
#' 
#' @examples
#' example_data <- system.file("extdata", "GSE117751_RAW", package = "NanoTube")
#' 
#' # Load data and positive control normalization
#' dat <- read_merge_rcc(list.files(example_data, full.names = TRUE))
#' dat <- normalize_pos_controls(dat)
#' 
#' # Remove endogenous genes that fail to reject the null hypothesis
#' # in a one-sided t test against negative control genes with p < 0.05.
#' dat <- remove_background(dat, mode = "t.test", pval = 0.05)
#' 
#' # Remove endogenous genes where fewer than 25% of samples have an expression
#' # 2 standard deviations above the average negative control gene. Also, 
#' # subtract this background level (mean + 2*sd) from endogenous genes.
#' dat <- remove_background(dat, mode = "threshold", 
#'                          numSD = 2, proportionReq = 0.25, subtract = TRUE)

remove_background <- function(dat, 
                              mode = c("threshold", "t.test"), 
                              numSD, proportionReq, pval,
                              subtract = FALSE) {
  
    # Check for negative control genes in the data set
    if (sum(dat$dict$CodeClass == "Negative") == 0) {
        stop("Cannot conduct filtering with negative controls: No negative 
            control genes found in input")
    }
    
    negative.mean <- apply(dat$exprs[dat$dict$CodeClass == "Negative",], 
                           2, mean)
    negative.sd <- apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, sd)
    
    if (mode == "threshold" | subtract) {
        bg.threshold <- negative.mean + numSD * negative.sd
        
        dat$bg.stats <- data.frame(Mean.Neg = negative.mean,
                                   Max.Neg =
                                     apply(dat$exprs[
                                       dat$dict$CodeClass == "Negative", ],
                                           2, max), 
                                   sd.Neg = negative.sd,
                                   background = bg.threshold,
                                   num.less.bg = 
                                     colSums(t(t(dat$exprs[
                                       dat$dict$CodeClass == "Endogenous" , ]) <
                                                 bg.threshold)), 
                                   frc.less.bg = 
                                     colMeans(t(t(dat$exprs[
                                       dat$dict$CodeClass == "Endogenous" , ]) <
                                                  bg.threshold)))
    } else {
        dat$bg.stats <- data.frame(Mean.Neg = negative.mean,
                                   Max.Neg = apply(dat$exprs[
                                     dat$dict$CodeClass == "Negative",], 
                                     2, max),
                                   sd.Neg = negative.sd)
    }
    
    if (mode == "t.test") {
        dat$gene.stats <- data.frame(row.names = dat$dict$Name,
                                     t.stat = rep(NA, 
                                                  times=length(dat$dict$Name)),
                                     p.val = rep(NA, 
                                                 times=length(dat$dict$Name)),
                                     pass = rep(NA, 
                                                times=length(dat$dict$Name)))
        
        rowsKeep <- ifelse(dat$dict$CodeClass == "Endogenous", 
                           yes = FALSE, no = TRUE)
        background <- unlist(dat$exprs[dat$dict$CodeClass == "Negative",])
        
        for (i in which(dat$dict$CodeClass == "Endogenous")){
            ttest <- t.test(x = dat$exprs[i,], y = background, 
                            var.equal = FALSE, alternative = "greater")
            if (ttest$p.value < pval) {
                rowsKeep[i] <- TRUE
            }
            dat$gene.stats$t.stat[i] <- as.numeric(ttest$statistic)
            dat$gene.stats$p.val[i] <- ttest$p.value
        }
        
        dat$gene.stats$pass <- rowsKeep
    } else {
        rowsKeep <- which(dat$dict$CodeClass != "Endogenous" | 
                            rowMeans(t(t(dat$exprs) > bg.threshold)) >= 
                            proportionReq)
    }
    
    if (subtract) {
        for (i in seq_len(ncol(dat$exprs))) dat$exprs[,i] <- 
            dat$exprs[,i] - bg.threshold[i]
        dat$exprs[dat$exprs < 0] <- 0
    }
    
    dat$exprs <- dat$exprs[rowsKeep,]
    dat$dict <- dat$dict[rowsKeep,]
    return(dat)
}
