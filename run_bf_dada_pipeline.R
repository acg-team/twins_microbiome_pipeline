#################### 0: CONFIGURATION : PLEASE SET!  ##############
conf <- vector(mode="list", length=3)
names(conf) <- c("location", "dataset", "pipeline")

### now we set it as calculare on local macbook and use 34 new dataset
conf$location <- "LOCAL"  # LOCAL / HOMESERVER  / ETHSERVER
conf$dataset <- "BODYFL"    #   TWIN / "BODYFL" /
conf$pipeline <- "DADA2"   # QIIME / DADA2

################# FILTERING parameters 
dada_param <- vector(mode="list", length=6)
names(dada_param) <- c("QUALITY_THRESHOLD", "maxEE", "trimLeft", "truncLen", "MSA_aligner", "tree_method")

dada_param$QUALITY_THRESHOLD <- 2
dada_param$maxEE <- c(2,2)
dada_param$trimLeft <- c(3,3)
dada_param$trimRight <- c(3,3)

dada_param$MSA_aligner <- "MUSCLE"   # DECIPHER  MUSCLE  clustalw 
dada_param$tree_method <- "RAXML"    # PHANGORN   

##################################################################

#TODO
# - add removal to QIIME2 and to TWIN pipeline


project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)


# full workflow (1-2 days on a server)
# cutadapt shall come here as well
#source("src/pipeline_dada2/1_metadata.R")

source("src/pipeline_dada2/2_file_names_parsing.R")

print("==================> long dada2 analysis has started...")
source("src/pipeline_dada2/4_BIG_dada_SV_table.R")

print("==================> Taxonomy assignment has started...")
source("src/pipeline_dada2/5_Tax_Assign.R")

print("==================> Phylogeny reconstraction has started...")
source("src/pipeline_dada2/6_Phylogeny.R")

source("src/pipeline_dada2/7_Create_Phyloseq_obj.R")

print(" Now PhyloSeq object has been created and you can run your analysis")
print(" >>>>>>>  END  <<<<<<<<")



