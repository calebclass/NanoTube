#' Unzip
#'
#' Unzips provided list of directories
#'
#' @param fileDirs character list of zip files
#' @return Unzipped directories

unzip.dirs <- function(fileDirs) {
  print(getwd())
  new.fileDirs <- gsub("\\.zip|\\.ZIP", "", fileDirs)
  for (i in 1:length(fileDirs)) {
    unzip(fileDirs[i], exdir = new.fileDirs[i])
  }
  return(new.fileDirs)
}
