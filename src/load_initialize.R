##### 0: load the necessary packages #####   

.cran_packages <- c("ggplot2", "gridExtra", "ShortRead", "XML", "tictoc", "MASS", "phangorn ")
.bioc_packages <- c("ShortRead", "devtools", "dada2", "phyloseq", "DECIPHER", "phangorn")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(suppressUpdates = FALSE)
  biocLite(.bioc_packages[!.inst], ask = F)
}

devtools::install_github("benjjneb/dada2")

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

set.seed(100)

#########################################################################

data_path    <- file.path(project_path, "/data/raw")
filt_path    <- file.path(data_path, "/data/processed")
result_path  <- file.path(project_path, "/reports")
rdp_path     <- file.path(project_path, "RDP")
silva_path   <- file.path(project_path, "/data/metadata")



