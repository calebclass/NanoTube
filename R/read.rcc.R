#' Read .rcc file
#'
#' This function reads in the .rcc files and splits into expression, sample data,
#' and qc components.
#'
#' @param file file name
#' @return list containing expression data, sample attributes, and basic qc from
#' the .rcc file.

read.rcc <- function(file) {

  dat.csv <- read.csv(file,
                      row.names = NULL, col.names = 1:4, header = FALSE, stringsAsFactors = FALSE)
  dat <- list(exprs = dat.csv[(which(dat.csv[,1] == "<Code_Summary>")+1):(which(dat.csv[,1] == "</Code_Summary>")-1),],
              sample = dat.csv[(which(dat.csv[,1] == "<Sample_Attributes>")+1):(which(dat.csv[,1] == "</Sample_Attributes>")-1), 1:2],
              qc = dat.csv[(which(dat.csv[,1] == "<Lane_Attributes>")+1):(which(dat.csv[,1] == "</Lane_Attributes>")-1), 1:2])
  for (i in 1:3) {
    colnames(dat[[i]]) <- dat[[i]][1,]
    dat[[i]] <- dat[[i]][-1,]
    rownames(dat[[i]]) <- NULL
  }
  dat[[1]]$Count <- as.numeric(dat[[1]]$Count)
  return(dat)
}
