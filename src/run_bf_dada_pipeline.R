################## run configuration ################################################
# it sets all paths accorting to parameters above
# TODO: use make instead?
source("src/load.R")
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)
source("src/configure.R")


#################### 0: CONFIGURATION : PLEASE SET!  ##############

### now we set it as calculare on local macbook and use 34 new dataset
conf$location <- "LOCAL"  # LOCAL / HOMESERVER  / ETHSERVER
conf$dataset <- "BFL"    #   TWIN / "BFL" /
conf$pipeline <- "DADA2"   # QIIME / DADA2

################# DADA FILTERING parameters 
dada_param$QUALITY_THRESHOLD <- 2
dada_param$maxEE <- c(2,4)

# trim If primers are at the start of your reads and are a constant length
dada_param$trimLeft <- c(0,0)
dada_param$trimRight <- c(0,0)

# be carefull, reads less then that are discarded!
dada_param$truncLen <-c(225,210)   # 230 / 210

################# MSA, Tree and Taxonomy parameters
tools_param$MSA_aligner <- "DECIPHER"   # DECIPHER  MUSCLE  clustalw 
tools_param$tree_method <- "RAXML"    # PHANGORN   
tools_param$tax_db <- "silva/silva_nr99_v138_train_set.fa.gz"  # "green_genes/gg_13_8_train_set_97.fa.gz"



########### Start pipeline
#TODO
# - add removal to QIIME2 and to TWIN pipeline



# cutadapt shall come here as well
#source("src/pipeline_dada2/1_metadata.R")

source("src/pipeline_dada2/2_file_names_parsing.R")

print("==================> long dada2 analysis has started...")
source("src/pipeline_dada2/4_BIG_dada_SV_table.R")

print("==================> Phylogeny reconstraction has started...")
source("src/pipeline_dada2/5_Phylogeny.R")


#source("src/pipeline_dada2/6_Tax_Assign_dada2_RDP.R")
source("src/pipeline_dada2/6_Tax_Assign_MapSeq.R")


print("==================> Creating the final results file...")
source("src/pipeline_dada2/7_Create_Phyloseq_obj.R")

print(" Now PhyloSeq object has been created and you can run your analysis")
print(" >>>>>>>  END  <<<<<<<<")



