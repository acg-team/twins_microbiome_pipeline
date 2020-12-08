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
source("src/load.R")
source("src/configure.R")
setwd(project_path)
load(file=file.path(files_intermediate_dada, seqtab.file)) 

# get a character vector of sequencess to assign taxonomy
seq <- dada2::getSequences(seqtab)
names(seq) <- seq
sequences <- Biostrings::DNAStringSet(seq)  #convert to DNAStringSet

# save sequences as a fasta file on disk, so MapSeq can use it
Biostrings::writeXStringSet(sequences, file=file.path(files_intermediate_dada, "sequences.fasta"))

# run MapSeq
# https://stackoverflow.com/questions/11395217/run-a-bash-script-from-an-r-script
# form a bash string
x1 <- 'mapseq -nthreads 4'
x2 <- file.path(files_intermediate_dada, "sequences.fasta")
x3 <- file.path(files_intermediate_dada, "sequences_taxonomy.fa.mseq")
x <- paste(x1,x2,'>',x3)
cat(x)

# run bash command
system(x)

# extract taxonomic classification from mseq file, convert to taxtab
mapseq.df <- read.table(x3, sep="\t", header=T, comment.char="", skip = 1)
rownames(mapseq.df) <- mapseq.df[,1]

taxtab <- mapseq.df[names(mapseq.df) %in% c("Phylum","Class","Order","Family","Genus")]

# save 
tax.fname <- paste0(substring(tools_param$tax_db, 1, 3), "_taxtab_mapseq.RData")


