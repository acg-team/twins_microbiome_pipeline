#### This script loads all neccesary packages and set up paths 
# Specific Packages used:
# ips -  Interfaces to Phylogenetic Software in R (to RAXML, Beast, MrBayes etc)
# dada2, phyloseq
# fastqcr - assessing FastQC quality in batch

# check environment
paste0(R.Version()[c("major","minor")], collapse = ".") # was 3.5.1
packageVersion("dada2") # 1.15
sessionInfo()

##### 0: load the necessary packages #####   
.cran_packages <- c("ggplot2", "gridExtra", "XML", "tictoc", "MASS", "ape", "phangorn",
                    "ips","tools","rphast","BiocManager","ggpubr", "tidyverse", "fastqcr", "plotly")

.bioc_packages <- c("ShortRead", "devtools", "dada2", "phyloseq", "msa", "DECIPHER")

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

#### 4: initialize configurarion data structure
conf <- vector(mode="list", length=3)
names(conf) <- c("location", "dataset", "pipeline")

dada_param <- vector(mode="list", length=5)
names(dada_param) <- c("QUALITY_THRESHOLD", "maxEE", "trimLeft", "trimRight", "truncLen")

tools_param <- vector(mode="list", length=4)
names(tools_param) <- c("MSA_aligner", "tree_method", "tax_db", "tax_method")


