#' Conduct differential expression analysis
#'
#' Use Limma to conduct a simple differential expression analysis. All groups 
#' are compared against the base.group, and empirical Bayes method is used to
#' identify significantly differentially expressed genes. Alternatively, a
#' design matrix can be supplied, as explained in limma::limmaUsersGuide()
#' 
#' @export
#'
#' @param dat NanoString data ExpressionSet, from processNanostringData
#' @param groups character vector, in same order as the samples in dat. NULL
#' if already included in 'dat'
#' @param base.group the group against which other groups are compared (must
#' be one of the levels in 'groups'). Will use the first group if NULL.
#' @param design a design matrix for Limma analysis (default NULL, will do 
#' analysis based on provided 'group' data)
#' @return The fit Limma object
#' 
#' @examples 
#' example_data <- system.file("extdata", "GSE117751_RAW", package = "NanoTube")
#' sample_info <- system.file("extdata", "GSE117751_sample_data.csv", 
#'                            package = "NanoTube")
#' 
#' dat <- processNanostringData(nsFiles = example_data,
#'                              sampleTab = sample_info, 
#'                              groupCol = "Sample_Diagnosis")
#' 
#' # Compare the two diseases against healthy controls ("None")
#' limmaResults <- runLimmaAnalysis(dat, base.group = "None")
#' 
#' 
#' # You can also supply a design matrix
#' # Generate fake batch labels
#' batch <- rep(c(0, 1), times = ncol(dat) / 2)
#' 
#' # Reorder groups ("None" first)
#' group <- factor(dat$groups, levels = c("None", "Autoimmune retinopathy", 
#'                                        "Retinitis pigmentosa"))
#' 
#' # Design matrix including sample group and batch
#' design <- model.matrix(~group + batch)
#' 
#' # Analyze data
#' limmaResults2 <- runLimmaAnalysis(dat, design = design)



runLimmaAnalysis <- function(dat, groups = NULL, base.group = NULL,
                             design = NULL) {

    dat.limma <- dat[fData(dat)$CodeClass == "Endogenous",]
    rownames(dat.limma) <- fData(dat)$Name[fData(dat)$CodeClass == "Endogenous"]
    
    # If RUV normalization was used, data are already log-transformed.
    if (dat$normalization[1] != "RUV") exprs(dat.limma) <- 
        log2(exprs(dat.limma) + 0.5)
  
    
    # Generate the design matrix for model fitting, if not provided
    if (is.null(design)) {
        # If group info was provided in function command, generate a 
        # phenoData table.
        if (!is.null(groups)) {
            sampData <- data.frame(row.names = colnames(dat.limma),
                                   groups = groups)
            
            phenoData(dat.limma) <- AnnotatedDataFrame(sampData)
        }
        
        # If group info was loaded from a sample table, it should be included 
        # in 'dat'.
        if (is.null(groups) & !("groups" %in% colnames(pData(dat)))) {
            stop("Groups not defined. Group identifiers must be loaded during
                 processNanostringData() or runLimmaAnalysis().")
        }
        
        # If base group was not defined, set it as the first group provided.
        if (is.null(base.group)) {
            base.group <- dat.limma$groups[1]
        } else {
            # Check validity of base.group
            if (!(base.group %in% dat.limma$groups)) {
                stop("'base.group' must be in 'groups'")
            }
        }
        
        if (is(dat.limma$groups, "factor")) {
            dat.limma$groups <- 
              factor(dat.limma$groups, 
                     levels = c(base.group, levels(dat.limma$groups)
                                [levels(dat.limma$groups) != base.group]))
        } else {
            dat.limma$groups <- 
              factor(dat.limma$groups, 
                     levels = c(base.group, unique(dat.limma$groups)
                                [unique(dat.limma$groups) != base.group]))
        }
        
        design <- model.matrix(~1 + dat.limma$groups)
        colnames(design) <- 
          gsub(" ", ".", gsub("dat.limma\\$groups", "", colnames(design)))
        colnames(design)[1] <- "Intercept"
      
    } else {
        # Check that design matrix is correct length
        if (nrow(design) != ncol(dat.limma)) {
            stop("Design matrix and expression data
                 do not have the same number of samples.")
        }
    }
    
    fit <- limma::lmFit(dat.limma, design)
    limmaFit <- limma::eBayes(fit)
    limmaFit$q.value <- limmaFit$p.value
    for (i in seq_len(ncol(limmaFit$q.value))) limmaFit$q.value[,i] <- 
        p.adjust(limmaFit$p.value[,i], method = "BH")
  
    limmaFit$eset <- dat.limma
    
    return(limmaFit)
}
