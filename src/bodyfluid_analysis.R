# @AlexY, created Jan 2020
# do analysis of Body Fluid dataset

#### init: load packages and set path
project_path <- "~/Projects_R/twins_microbiome_pipeline"
files_intermediate    <- file.path(project_path, "data_set_bodyfl/files_intermediate")
source("src/load.R")
setwd(project_path)

### LOAD Budy Fluid PhyloSeq file + metadata
metadata_path   <- file.path(project_path, "data_set_bodyfl/metadata")
load(file=file.path(metadata_path, "metadata.RData"))
load(file=file.path(files_intermediate, "phyloseq_object.RData"))

# sanity check of loaded data
get_taxa_unique(ps.tweens, "Phylum")
get_taxa_unique(ps.tweens, "Genus")
# Q: is NA a problem?


