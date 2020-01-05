# datasets are switched here

#################### 0: CONFIGURATION : PLEASE SET!  ##############
conf <- vector(mode="list", length=2)
names(conf) <- c("location", "dataset")

### now we set it as calculare on local macbook and use 34 new dataset
conf$location <- "LOCAL"  # LOCAL / HOMESERVER  / ETHSERVER
conf$dataset <- "BODYFL"  #  TWIN / "BODYFL" /
##################################################################




###### 1: Set Paths to folders and files ###############################

# TODO create if they dont exists
#if(!file_test("-d", filt_path)) dir.create(result_path)
#if(!file_test("-d", filt_path)) dir.create(filt_path)  # create filetered folder

# set project path depending on location
if(conf$location == "LOCAL"){
  project_path <- "~/Projects_R/twins_microbiome_pipeline"
} else if(conf$location == "HOMESERVER") {
  project_path <- "/media/alex/db5547c3-1ac1-4ec5-aac9-29a383a87978"
} else {
  stop(" WRONG SERVER CONFIGURATION")
}


######## 2:  set pathes to folders depending on dataset
if(conf$dataset == "TWIN"){
  metadata_path   <- file.path(project_path, "data_set_twin/metadata")
  files_intermediate    <- file.path(project_path, "data_set_twin/files_intermediate")
  result_path  <- file.path(project_path, "data_set_twin/reports_generated")
  if(conf$location == "LOCAL"){
    data_path <- file.path(project_path, "data/raw")
    filt_path <- file.path(project_path, "data/raw/filtered")
  } else if(conf$location == "HOMESERVER"){
    data_path    <- file.path(project_path, "BIOINF_DATA/TwinUK_Full")
    filt_path    <- file.path(project_path, "BIOINF_DATA/TwinUK_Full/filtered")
  }
} else if (conf$dataset == "BODYFL"){
  metadata_path   <- file.path(project_path, "data_set_bodyfl/metadata")
  files_intermediate    <- file.path(project_path, "data_set_bodyfl/files_intermediate")
  result_path  <- file.path(project_path, "data_set_bodyfl/reports_generated")
  if(conf$location == "LOCAL"){
    data_path <- file.path(project_path, "data_set_bodyfl/raw")
    filt_path <- file.path(project_path, "data_set_bodyfl/raw/filtered")
  } else if(conf$location == "HOMESERVER"){
    stop("no such configuration exists")
  }
}



##### 3: Set file names (same for any dataset but in different folders) ###########################

# file names for intermediate results
metadata.file <- "metadata.RData"
dada.err.file <- "dada_err_data.RData"
mergers.file <- "mergers.RData"
seqtab.file <- "seqtab.RData"
seqtab.snames.file <- "seqtab_snames.RData"

taxtab.file <- "taxtab.RData"

msa.file <- "msa.RData"
phylo.file <- "phylo_trees.RData"

phyloseq.file <- "phyloseq_object.RData"
phyloseq_analysis.file <- "phyloseq_analysis.RData"
