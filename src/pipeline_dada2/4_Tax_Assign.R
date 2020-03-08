# Author: @Alex
# Project: UK Twins
# Taxonomy assignment to Seq Table
# http://benjjneb.github.io/dada2/assign.html
# http://benjjneb.github.io/dada2/tutorial.html

# TODO: it is told that silva assignment is not very accurate, please check


#### init: load packages and set path
source("src/load.R")
source("src/configure.R")
setwd(project_path)

#########################################################################

### LOAD PREVIOUS DATA
load(file=file.path(files_intermediate, mergers.file)) 
load(file=file.path(files_intermediate, seqtab.file)) 


###### 1: ASSIGN TAXONOMY #################################
# ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
###########################################################
# assignTaxonomy implements the RDP Naive Bayesian Classifier algorithm
# taxtab[1:1227, 1:6] - for each inferred sequence we have a taxonomy assignment

tic()
#ref_fasta <- file.path(rdp_path, "rdp_train_set_14.fa.gz")

# TODO what is the dufference btw training and assignment?
ref_fasta <- file.path(silva_path, "silva_nr_v132_train_set.fa.gz")
#ref_fasta <- file.path(metadata_path, "silva_species_assignment_v128.fa.gz")
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

save(taxtab, file=file.path(files_intermediate, taxtab.file)) 



###### 1: ASSIGN SPECIES #################################
# assignSpecies, is intended for the more difficult task of taxonomic assignments 
# at or near the resolution limit of the sequenced marker-gene
# Do we need to assign species?



