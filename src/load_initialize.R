# Packages used
# ips -  https://www.rdocumentation.org/packages/ips/versions/0.0.11

##### 0: load the necessary packages #####   
# ips - Interfaces to Phylogenetic Software in R (to RAXML, Beast, MrBayes etc)
.cran_packages <- c("ggplot2", "gridExtra", "ShortRead", "XML", "tictoc", "MASS", "ape", "phangorn","ips","tools")
.bioc_packages <- c("ShortRead", "devtools", "dada2", "phyloseq", "msa")

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

#devtools::install_github("benjjneb/dada2")

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)


# Check versions if nessesary
packageVersion("dada2")
R.version
sessionInfo()

set.seed(100)


###### 1: Paths to folders ###############################

# TODO create if they dont exists
#if(!file_test("-d", filt_path)) dir.create(result_path)
#if(!file_test("-d", filt_path)) dir.create(filt_path)  # create filetered folder

data_path    <- file.path(project_path, "data/raw")
filt_path    <- file.path(project_path, "data/raw/filtered")
processed_path <- file.path(project_path, "data/processed")
files_intermediate    <- file.path(project_path, "files_intermediate")

result_path  <- file.path(project_path, "reports_generated")
rdp_path     <- file.path(project_path, "RDP")
metadata_path   <- file.path(project_path, "data/metadata")


##### 2: Set parameters ###########################$#
# TODO: set trimming parameters in a single dataset (tuple?)

# file names for intermediate results
metadata.file <- "metadata.RData"
dada.err.file <- "dada_err_data_q15.RData"
mergers.file <- "mergers_q15.RData"
seqtab.file <- "seqtab_q15.RData"
seqtab.snames.file <- "seqtab_snames_q15.RData"
taxtab.file <- "taxtab_g15.RData"

msa.file <- "msa.RData"
phylo.file <- "fitGTR.RData"
treeGTR_2.file <- "fitGTR_2.RData"

phyloseq.file <- "phyloseq_g15.RData"
