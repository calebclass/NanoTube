#' Conduct differential expression analysis
#'
#' Use Limma to conduct a simple differential expression analysis. All groups are
#' compared against the base.group, and empirical Bayes method is used to
#' identify significantly differentially expressed genes.
#'
#' @param dat NanoString data, including expression matrix and gene dictionary
#' @param groups character vector, in same order as the samples in dat
#' @param base.group the group against which other groups are compared (must
#' be one of the levels in 'groups')
#' @return The fit Limma object

runLimmaAnalysis <- function(dat, groups, base.group) {

  if (!(base.group %in% groups)) stop("'base.group' must be in 'groups'")

  dat.limma <- dat$exprs[dat$dict$CodeClass == "Endogenous",]
  rownames(dat.limma) <- dat$dict$Name[dat$dict$CodeClass == "Endogenous"]
  dat.limma <- log2(dat.limma + 0.5)
  #dat.limma <- dat.limma[,order(colnames(dat.limma))]

  if (class(groups) == "factor") {
    groups.f <- factor(groups, levels = c(base.group, levels(groups)[levels(groups) != base.group]))
  } else {
    groups.f <- factor(groups, levels = c(base.group, unique(groups)[unique(groups) != base.group]))
  }
  sampData <- data.frame(row.names = colnames(dat.limma),
                         group = groups.f)

  exprSet <- ExpressionSet(assayData = as.matrix(dat.limma),
                           phenoData = AnnotatedDataFrame(sampData))

  design <- model.matrix(~0 + exprSet$group)
  colnames(design) <- gsub(" ", ".", gsub("exprSet\\$group", "", colnames(design)))

  #Set vehicle to intercept
  design[,1] <- 1
  colnames(design)[1] <- "Intercept"

  fit <- lmFit(exprSet, design)
  limmaFit <- eBayes(fit)
  limmaFit$q.value <- limmaFit$p.value
  for (i in 1:ncol(limmaFit$q.value)) limmaFit$q.value[,i] <- p.adjust(limmaFit$p.value[,i], method = "BH")

  limmaFit$eset <- exprSet
  limmaFit$dict <- dat$dict[dat$dict$CodeClass == "Endogenous",]
  limmaFit$sampleData <- sampData

  return(limmaFit)
}
