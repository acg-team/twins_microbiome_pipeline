# Taxonomy assignment with MapSeq classifier
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

#### init: load packages and set path
load(file=file.path(files_intermediate_dada, seqtab.file)) 

# get a character vector of sequencess to assign taxonomy
seq <- asv_sequences #dada2::getSequences(seqtab)

sequences.dna.stringset <- Biostrings::DNAStringSet(seq)  #convert to DNAStringSet

# save sequences as a fasta file on disk, so MapSeq can use it
Biostrings::writeXStringSet(sequences.dna.stringset, file=file.path(files_intermediate_dada, "asv_sequences.fasta"))


# run MapSeq
# https://stackoverflow.com/questions/11395217/run-a-bash-script-from-an-r-script
# form a bash string
x1 <- 'mapseq -nthreads 5 -minscore 0'
x2 <- file.path(files_intermediate_dada, "asv_sequences.fasta")
x3 <- file.path(files_intermediate_dada, "sequences_taxonomy.fa.mseq")
x <- paste(x1,x2,'>',x3)
cat(x)

# run bash command
system(x)

# extract taxonomic classification from mseq file, convert to taxtab (skip the first line with text)
mapseq.df <- read.table(x3, sep="\t", header=T, comment.char="", skip = 1)
rownames(mapseq.df) <- mapseq.df[,1]
#mapseq.df <- mapseq.df[sort(as.numeric(mapseq.df$X.query)),]

taxtab <- mapseq.df[names(mapseq.df) %in% c("Phylum","Class","Order","Family","Genus")]
#TODO : do we need to filer by score - id bitscore < 20, set NA


# SANITY CHECK: here we need to compare rownames(mapseq.df) and names(seq), should be zero difference
setdiff(names(seq), rownames(mapseq.df))



# save 
tax.fname <- paste0("taxtab_", tools_param$tax_db, "_mapseq.RData")
fname = file.path(files_intermediate_dada, tax.fname) 
save(taxtab, file=fname)

print(paste("Taxonomy saved to ", fname))

