#' Unzip
#'
#' Unzips provided list of directories
#' 
#' @export
#'
#' @param fileDirs character list of zip files
#' @return Names of now-unzipped directories

unzip_dirs <- function(fileDirs) {
    print(getwd())
    new.fileDirs <- gsub("\\.zip|\\.ZIP", "", fileDirs)
    for (i in seq_along(fileDirs)) {
        unzip(fileDirs[i], exdir = new.fileDirs[i])
    }
    return(new.fileDirs)
}
