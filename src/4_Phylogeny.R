#### init: load packages and set path
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)
source("src/load_initialize.R")

load(file=file.path(models_path, seqtab.file)) 
load(file=file.path(models_path, seqtab.snames.file)) 
load(file=file.path(models_path, taxtab.file))


######## Simple Phylogeny (MSA then NJ)  #######
#####  MSA 
seqs <- dada2::getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree

alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor=NA)

#### Construct tree 
### NOTE: if we use seqtab, this is not a tree of species but a tree of sequence variants!

tic()
# neibor join tree
phang.align <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
fit = phangorn::pml(treeNJ, data=phang.align)

######  GTR tree
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
toc()

save(fitGTR,treeNJ, file=file.path(models_path, "fitGTR.RData")) 

