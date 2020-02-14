install.packages(c("devtools", "plotly", "mHG"))

devtools::install_github("r-lib/rlang")
devtools::install_github("tidyverse/ggplot2")

source("https://bioconductor.org/biocLite.R")
biocLite("Biobase")
biocLite("limma")

devtools::install("NanoTube") #local install