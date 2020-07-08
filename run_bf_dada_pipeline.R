#################### 0: CONFIGURATION : PLEASE SET!  ##############
conf <- vector(mode="list", length=3)
names(conf) <- c("location", "dataset", "pipeline")

### now we set it as calculare on local macbook and use 34 new dataset
conf$location <- "LOCAL"  # LOCAL / HOMESERVER  / ETHSERVER
conf$dataset <- "BFL"    #   TWIN / "BFL" /
conf$pipeline <- "DADA2"   # QIIME / DADA2

################# FILTERING parameters 
dada_param <- vector(mode="list", length=5)
names(dada_param) <- c("QUALITY_THRESHOLD", "maxEE", "trimLeft", "MSA_aligner", "tree_method")

dada_param$QUALITY_THRESHOLD <- 2
dada_param$maxEE <- c(4,5)
dada_param$trimLeft <- c(3,3)
dada_param$trimRight <- c(3,5)

dada_param$MSA_aligner <- "DECIPHER"   # DECIPHER  MUSCLE  clustalw 
dada_param$tree_method <- "RAXML"    # PHANGORN   

##################################################################

#TODO
# - add removal to QIIME2 and to TWIN pipeline

project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)


# cutadapt shall come here as well
#source("src/pipeline_dada2/1_metadata.R")

source("src/pipeline_dada2/2_file_names_parsing.R")

print("==================> long dada2 analysis has started...")
source("src/pipeline_dada2/4_BIG_dada_SV_table.R")



print("==================> Taxonomy assignment has started...")
source("src/pipeline_dada2/5_Tax_Assign.R")

print("==================> Phylogeny reconstraction has started...")
source("src/pipeline_dada2/6_Phylogeny.R")

print("==================> Creating the final results file...")
source("src/pipeline_dada2/7_Create_Phyloseq_obj.R")

print(" Now PhyloSeq object has been created and you can run your analysis")
print(" >>>>>>>  END  <<<<<<<<")



