###### Taxonomy assignment with RDP Naive Bayesian Classifier built-in to dada2
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

print("==================> Taxonomy assignment has been started...")

#### init: load packages and set path
load(file=file.path(files_intermediate_dada, seqtab.file)) 

# get a character vector of sequencess to assign taxonomy
sequences <- asv_sequences #dada2::getSequences(seqtab)

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
cat("=> Total time of taxonomy assignment: ")
toc() # 1468 sec = 24 min


# modify the output format for Green Genes
# Green Genes adds a 7th column (species) and also aff f_ in front of all taxa, remove it
if (tools_param$tax_db == "green_genes/gg_13_8_train_set_97.fa.gz"){
  taxtab <- taxtab[,-7]
  taxtab <- apply(taxtab, 1:2, gsub, pattern='^(k|p|c|o|f|g)__', replacement='')
}



# Save to disk
# TODO: do it as csv? [sequence, phylim, kingdom ...]
tax.fname <- paste0("taxtab_", substring(tools_param$tax_db, 1, 3), "_dadardp.RData")
fname = file.path(files_intermediate_dada, tax.fname) 
save(taxtab, file=fname)
print(paste("saved to ", fname))





