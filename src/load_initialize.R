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

# TODO create if they dont exists
#if(!file_test("-d", filt_path)) dir.create(result_path)
#if(!file_test("-d", filt_path)) dir.create(filt_path)  # create filetered folder

data_path    <- file.path(project_path, "data/raw")
filt_path <- file.path(project_path, "data/raw/filtered")
processed_path    <- file.path(project_path, "data/processed")
models_path    <- file.path(project_path, "models")

result_path  <- file.path(project_path, "reports")
rdp_path     <- file.path(project_path, "RDP")
metadata_path   <- file.path(project_path, "data/metadata")



