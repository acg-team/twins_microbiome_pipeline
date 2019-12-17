# Packages used
# ips -  https://www.rdocumentation.org/packages/ips/versions/0.0.11


# TODO change to BiocManager /  https://www.bioconductor.org/packages/release/bioc/html/dada2.html

paste0(R.Version()[c("major","minor")], collapse = ".") # was 3.5.1
packageVersion("dada2") # 1.15
sessionInfo()

##### 0: load the necessary packages #####   
# ips - Interfaces to Phylogenetic Software in R (to RAXML, Beast, MrBayes etc)
.cran_packages <- c("ggplot2", "gridExtra", "XML", "tictoc", "MASS", "ape", "phangorn","ips","tools","rphast","BiocManager")
.bioc_packages <- c("ShortRead", "devtools", "dada2", "phyloseq", "msa")

# CRAN packages
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
} else{
  print(" All CRAN packages are already installed.")
}

# Bioconductor packages
BiocManager::version()
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  BiocManager::install(.bioc_packages[!.inst])
} else{
  print(" All Bioconductor packages are already installed.")
  # use BiocManager::install() if you need to update all
}

#devtools::install_github("benjjneb/dada2")

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)



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
phyloseq_analysis.file <- "phyloseq_analysis.RData"
