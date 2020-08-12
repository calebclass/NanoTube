#' Convert NanoString ExpressionSet to NanoStringSet
#'
#' Convert ExpressionSet from processNanoStringData to a NanoStringSet for use
#' in NanoStringDiff.
#'
#' @param eset NanoString data ExpressionSet, from processNanostringData
#' @param designs Design matrix. If NULL, will look for "groups" column in 
#' pData(eset).
#' @return A NanoStringSet for NanoStringDiff

createNanoStringSetFromEset <- function(eset, designs = NULL) {
  
  # Order positive control data (required for NanoStringDiff)
  pos.unordered <- eset[fData(eset)$CodeClass == "Positive",]
  pos.ordered <- exprs(pos.unordered)[order(fData(pos.unordered)$Name),]
  
  # Look for group info
  if (is.null(designs)) {
    if ("groups" %in% colnames(pData(eset))) {
      designs <- model.matrix(~0 + eset$groups)
    } else {
      stop("Must input design matrix, or include groups in sample info of eset.")
    }
  }
  
  # Build NanoStringSet
  nsSet <- NanoStringDiff::createNanoStringSet(endogenous = exprs(eset)[fData(eset)$CodeClass == "Endogenous",],
                                               positiveControl = pos.ordered,
                                               negativeControl = exprs(eset)[fData(eset)$CodeClass == "Negative",],
                                               housekeepingControl = exprs(eset)[fData(eset)$CodeClass == "Housekeeping",],
                                               designs = designs)
  
  return(nsSet)
}