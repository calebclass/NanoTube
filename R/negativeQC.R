negativeQC <- function(nsRaw, numSD = 2) {
  
  if (nsRaw$normalization[1] != "none") {
    stop("Raw NanoString data must be provided. Use normalization='none'
         in `processNanostringData`.")
  }
  
  dat.neg <- as.data.frame(nsRaw$dat$exprs.raw[nsRaw$dat$dict.raw$CodeClass == "Negative",])
  
  # Statistics for negative control genes
  negative.mean <- apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, mean)
  negative.sd <- apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, sd)
  bg.threshold <- negative.mean + numSD * negative.sd
  
  # Summary Table
  dat$bg.stats <- data.frame(Mean.Neg = negative.mean,
                             Max.Neg = apply(dat$exprs[dat$dict$CodeClass == "Negative",], 2, max),
                             sd.Neg = negative.sd,
                             background = bg.threshold,
                             num.less.bg = colSums(t(t(dat$exprs[dat$dict$CodeClass == "Endogenous" ,]) < bg.threshold)),
                             frc.less.bg = colMeans(t(t(dat$exprs[dat$dict$CodeClass == "Endogenous" ,]) < bg.threshold)))
  neg.tab$fail <- paste0(nsRaw$bg.stats$num.less.bg, " (",
                         round(nsRaw$bg.stats$frc.less.bg*100, 1), "%)")
  
  neg.tab <- neg.tab[,-(5:6)]
  colnames(neg.tab) <- c("Mean (Neg)", "Max (Neg)", "sd (Neg)", "BG (Mean+nSD)", 
                         "Genes below BG (%)")
  
  
  # Strip plot for negative control genes
  dat.neg$Gene <- nsRaw$dat$dict.raw$Name[nsRaw$dat$dict.raw$CodeClass == "Negative"]
  
  dat.neg.df <- reshape::melt(dat.neg, "Gene")
  colnames(dat.neg.df)[2:3] <- c("Sample", "Count")

  
  neg1 <- ggplot(data = dat.neg.df, aes(x=Count, y=Sample, 
                                        text=paste0("Sample: ", Sample, "\nGene: ", Gene, "\nCount: ", Count))) +
    geom_jitter(height = 0.2, width = 0, colour = "black", fill = "grey70", pch=21) +
    theme_classic() + ylab("") 
  
  #neg.plot <- div(ggplotly(neg1, tooltip = c("text"), width = 550, height = 400) %>% 
  #            layout(margin = list(l=90), autosize = FALSE), 
  #                    align="center")
  neg.plot <- ggplotly(neg1, tooltip = c("text"), width = 550, height = 400) %>% 
    layout(margin = list(l=90), autosize = FALSE)
  
  return(list(tab = neg.tab,
              plt = neg.plot))
}