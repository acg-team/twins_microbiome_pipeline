################## run configuration 
# it sets all paths accorting to parameters above
# TODO: use make instead?
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)

# load packages and initialize config variables
source("src/load.R")


# configure paths depending on server / datasets etc 
conf$location <- "LOCAL"  # LOCAL / HOMESERVER  / ETHSERVER
conf$dataset <- "BFL"    #   TWIN / "BFL" /
conf$pipeline <- "DADA2"   # QIIME / DADA2
source("src/configure.R")


########### Start pipeline #############
# - add removal to QIIME2 and to TWIN pipeline


# cutadapt shall come here as well
#source("src/pipeline_dada2/1_metadata.R")

source("src/pipeline_dada2/2_file_names_parsing.R")


print("==================> long dada2 analysis has started...")
dada_param$QUALITY_THRESHOLD <- 2
dada_param$maxEE <- c(2,4)

# trim If primers are at the start of your reads and are a constant length
dada_param$trimLeft <- c(0,0)
dada_param$trimRight <- c(0,0)

# be carefull, reads less then that are discarded!
dada_param$truncLen <-c(210,220)   # 230 / 210
# INPUT:
source("src/pipeline_dada2/4_BIG_dada_SV_table.R")
# OUTPUT:


print("==================> Phylogeny reconstraction has started...")
# INPUT:
tools_param$MSA_aligner <- "DECIPHER"   # DECIPHER / MUSCLE / clustalw 
tools_param$tree_method <- "RAXML"    # PHANGORN   
source("src/pipeline_dada2/5_Phylogeny.R")
# OUTPUT:


print("==================> Taxonomy assignment has started...")
# INPUT:
#tools_param$tax_db <- "silva/silva_nr99_v138_train_set.fa.gz"  # "green_genes/gg_13_8_train_set_97.fa.gz"
tools_param$tax_db <- "ncbi"
tools_param$tax_method <- "mapseq"
#source("src/pipeline_dada2/6_Tax_Assign_dada2_RDP.R")
source("src/pipeline_dada2/6_Tax_Assign_MapSeq.R")
# OUTPUT


print("==================> Creating the final results file...")
source("src/pipeline_dada2/7_Create_Phyloseq_obj.R")


print(" Now PhyloSeq object has been created and you can run your analysis")
print(" >>>>>>>  END  <<<<<<<<")



