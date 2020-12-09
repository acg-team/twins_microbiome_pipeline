print(paste("########   Configuration set for : ", conf, " ###############" ))


###### 1: Set Paths to folders and files ###############################
# TODO create if they dont exists
#if(!file_test("-d", filt_path)) dir.create(result_path)
#if(!file_test("-d", filt_path)) dir.create(filt_path)  # create filetered folder


# General Project and installed utilities 
if(conf$location == "LOCAL"){
  project_path <- "~/Projects_R/twins_microbiome_pipeline"
  raxm.exec.path <- "/Users/alex/bioinf_tools/RAxML/raxmlHPC-PTHREADS-AVX"
  fastqc.path = "/Applications/BIOINF/FastQC.app/Contents/MacOS/fastqc"

} else if(conf$location == "HOMESERVER") {
  project_path <- "/media/alex/db5547c3-1ac1-4ec5-aac9-29a383a87978"
  raxm.exec.path <- "/home/alex/installed/BIOINF_tools/RAxML/raxmlHPC-PTHREADS-AVX"
  
} else if(conf$location == "NATASHA") {
  project_path <- "/Users/Natasha/Mac_Documents/Academic_Research/Forensics/Projects/Twin_Microbiome_Data/twins_microbiome_pipeline-master"
  raxm.exec.path <- "????/home/alex/installed/BIOINF_tools/RAxML/raxmlHPC-PTHREADS-AVX"

} else {
  stop(" WRONG SERVER CONFIGURATION")
}


taxonomy_db_path <- file.path(project_path, "16S_taxonomy_db")


######## 2:  set SPECIFIC FOLDERS depending on dataset
if(conf$dataset == "TWIN"){
  metadata_path   <- file.path(project_path, "data_set_twin/metadata")
  files_intermediate_dada  <- file.path(project_path, "data_set_twin/files_intermediate_dada")
  files_intermediate_qiime <- file.path(project_path, "data_set_twin/files_intermediate_qiime")
  result_path  <- file.path(project_path, "data_set_twin/reports_generated")
  
  if(conf$location == "LOCAL"){
    data_path <- file.path(project_path, "data_set_twin/fastq")
    filt_path <- file.path(project_path, "data_set_twin/fastq/filtered")
    qiime_qza_path <- file.path(project_path, "data_set_twin/fastq/qza")
  } else if(conf$location == "HOMESERVER"){
    data_path    <- file.path(project_path, "BIOINF_DATA/TwinUK_Full")
    filt_path    <- file.path(project_path, "BIOINF_DATA/TwinUK_Full/filtered")
  }

  
} else if (conf$dataset == "BFL"){
  metadata_path   <- file.path(project_path, "data_set_bodyfl/metadata")
  files_intermediate_dada  <- file.path(project_path, "data_set_bodyfl/files_intermediate_dada")
  files_intermediate_qiime <- file.path(project_path, "data_set_bodyfl/files_intermediate_qiime")
  result_path  <- file.path(project_path, "data_set_bodyfl/reports_generated")

  if(conf$location == "LOCAL"){
    raw_data_path <- file.path(project_path, "data_set_bodyfl/fastq/no_primers")
    filt_path <- file.path(project_path, "data_set_bodyfl/fastq/filtered")
    qiime_qza_path <- file.path(project_path, "data_set_bodyfl/fastq/qza")
  } else if(conf$location == "HOMESERVER"){
    stop("no such configuration exists yet")
  }
}



##### 3: Set file names (same for any dataset but in different folders) ###########################

### configuration shall be set before!
# form file names for intermediate results
metadata.file <- "metadata.RData"
dada.err.file <- "dada_err_data.RData"
mergers.file <- "mergers.RData"
seqtab.file <- "seqtab.RData"   #paste0("seqtab_", file.suffix)
seqtab.snames.file <- "seqtab_snames.RData"  #paste0("seqtab_snames_", file.suffix)

msa.file <- "msa.RData"
phylo.file <- "phylo_trees.RData"





