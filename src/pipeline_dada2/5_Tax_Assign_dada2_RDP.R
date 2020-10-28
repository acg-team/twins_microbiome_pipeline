# Taxonomy assignment with RDP Naive Bayesian Classifier built-in to dada2
# http://benjjneb.github.io/dada2/assign.html
# http://benjjneb.github.io/dada2/tutorial.html
# Author: @AlexY

# INPUT: 
#  - one of the reference 16S databases: silva, green genes
#  - a vector of sequences to classify (extract from seqtab)

# OUTPUT
# - taxtab, a matrix of sequences x ["Kingdom", "Phylum" ... ]


###### ASSIGN TAXONOMY ####################################
# ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
###########################################################

print("==================> Taxonomy assignment has started...")

#### init: load packages and set path
source("src/load.R")
source("src/configure.R")
setwd(project_path)
load(file=file.path(files_intermediate_dada, seqtab.file)) 

# get a character vector of sequencess to assign taxonomy
sequences <- dada2::getSequences(seqtab)

### Choose a reverence database (SILVA, green genes, RDP)
ref_fasta <- file.path(taxonomy_db_path, tools_param$tax_db)
print(paste("Start assignment with taxonomy db : ", tools_param$tax_db))


# assignTaxonomy implements the RDP Naive Bayesian Classifier algorithm
# taxtab[1:1227, 1:6] - for each inferred sequence we have a taxonomy assignment
tic()
taxtab <- dada2::assignTaxonomy(
  seqs = sequences, 
  refFasta = ref_fasta,
  multithread = 3,
  minBoot=50
  )
print("Total time of taxonomy assignment:")
toc() # 1468 sec = 24 min


# modify the output format for Green Genes
# Green Genes adds a 7th column (species) and also aff f_ in front of all taxa, remove it
if (tools_param$tax_db_name){
  taxtab <- taxtab[,-7]
  taxtab <- apply(taxtab, 1:2, gsub, pattern='^(k|p|c|o|f|g)__', replacement='')
}



# Save to disk
# TODO: do it as csv? [sequence, phylim, kingdom ...]
fname = file.path(files_intermediate_dada, paste0(tools_param$tax_db_name, "_taxtab.RData")) 
save(taxtab, file=fname)
print(paste("saved to ", fname))






######  3 - Kraken





