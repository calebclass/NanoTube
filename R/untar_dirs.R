#' Untar
#'
#' Untars provided list of directories (analogous to unzip_dirs)
#' 
#' @export
#'
#' @param fileDirs character list of tar files
#' @return Names of now-untarred directories

untar_dirs <- function(fileDirs) {
    print(getwd())
    new.fileDirs <- gsub("\\.tar|\\.TAR", "", fileDirs)
    for (i in seq_along(fileDirs)) {
        untar(fileDirs[i], exdir = new.fileDirs[i])
    }
    return(new.fileDirs)
}
