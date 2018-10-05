# @Alex: UK Twins
# Taxonomy assignment to Seq Table
# http://benjjneb.github.io/dada2/assign.html
# http://benjjneb.github.io/dada2/tutorial.html


#### init: load packages and set path
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)
source("src/load_initialize.R")

#########################################################################

### LOAD PREVIOUS DATA
load(file=file.path(models_path, mergers.file)) 
load(file=file.path(models_path, seqtab.file)) 


###### 1: ASSIGN TAXONOMY ########################
# assignTaxonomy implements the RDP Naive Bayesian Classifier algorithm
# taxtab[1:1227, 1:6] - for each inferred sequence we have a taxonomy assognment

tic()
#ref_fasta <- file.path(rdp_path, "rdp_train_set_14.fa.gz")
ref_fasta <- file.path(metadata_path, "silva_nr_v128_train_set.fa.gz")
taxtab <- dada2::assignTaxonomy(seqtab, refFasta = ref_fasta)
colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
unname(taxtab)
print("Total time of taxonomy assignment:")
toc()

## now we have
# - seqtab [1:28, 1:1227]
# - taxtab [1:1227, 1:6]
# so we can combine for each sample a table with taxa names and abanduncies
# see Exploratoty_Analysys file

save(taxtab, file=file.path(result_path, taxtab.file)) 

