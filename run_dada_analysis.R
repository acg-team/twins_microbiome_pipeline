#################### 0: CONFIGURATION : PLEASE SET!  ##############
conf <- vector(mode="list", length=3)
names(conf) <- c("location", "dataset", "pipeline")

### now we set it as calculare on local macbook and use 34 new dataset
conf$location <- "LOCAL"  # LOCAL / HOMESERVER  / ETHSERVER
conf$dataset <- "BODYFL"    #   TWIN / "BODYFL" /
conf$pipeline <- "DADA2"   # QIIME / DADA2

################# FILTERING parameters 
dada_param <- vector(mode="list", length=2)
names(dada_param) <- c("QUALITY_THRESHOLD", "maxEE")
dada_param$QUALITY_THRESHOLD <- 2
dada_param$maxEE <- c(2,4)
dada_param$MSA_aligner <- "MUSCLE"   # DECIPHER  MUSCLE  clustalw 
dada_param$tree_method <- "RAXML"    # PHANGORN   
##################################################################



#TODO
# - remove primers and adapters - done
# - add removal to QIIME2 
# - keep track of reads filtered before and after - done
# - add exception in case merging is failed

project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)

#con <- file("last_dada.log")
#sink(con, append=TRUE)
#sink(con, append=TRUE, type="message")

# full workflow (1-2 days ona  a server)
#source("src/pipeline_dada2/1_metadata.R")
source("src/pipeline_dada2/2_file_names_parsing.R")

print("==================> long dada2 analysis has started...")
source("src/pipeline_dada2/3_BIG_dada_SV_table.R")

print("==================> Taxonomy assignment has started...")
source("src/pipeline_dada2/4_Tax_Assign.R")

print("==================> Phylogeny reconstraction has started...")
source("src/pipeline_dada2/5_Phylogeny.R")

source("src/pipeline_dada2/6_Create_Phyloseq_obj.R")

print(" Now PhyloSeq object has been created and you can run your analysis")
print(" >>>>>>>  END  <<<<<<<<")
# Restore output to console
#sink() 
#sink(type="message")


