#### INFER A PHYLOGENY TREE ( beta diversity?)
##############################################

#### init: load packages and set path
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)
source("src/load_initialize.R")

load(file=file.path(files_intermediate, seqtab.file)) 
load(file=file.path(files_intermediate, seqtab.snames.file)) 
load(file=file.path(files_intermediate, taxtab.file))


######## Simple Phylogeny (MSA then NJ)  #######
#####  MSA 
seqs <- dada2::getSequences(seqtab)
names(seqs) <- seqtab.samples.names # This propagates to the tip labels of the tree

#msa package provides a unified R/Bioconductor interface to MSA: ClustalW, ClustalOmega, and Muscle
# TODO: check for ClustalW specific parameters
microbiome.msa.clustalW <- msa::msa(seqs, method="ClustalW", type="dna", order="input")

# TODO: check for Muscle specific parameters
microbiome.msa.muscle <- msa::msa(seqs, method="Muscle", type="dna", order="input")

# save MSA to a file 
save(microbiome.msa.clustalW,microbiome.msa.muscle, file=file.path(files_intermediate, msa.file)) 




########## Construct tree 
### NOTE: if we use seqtab, this is not a tree of species but a tree of sequence variants!

## Option 1: GTR tree
tic()
# infer a guide tree
phang.align <- as.phangorn::phyDat(mult, type="DNA", names=seqtab.samples.names)
dm <- dist.ml(phang.align)
treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
fit = phangorn::pml(treeNJ, data=phang.align)

## infer a  GTR tree
fitGTR <- update(fit, k=4, inv=0.2)  #???
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
toc()
#########################


# save the tree to file
save(fitGTR, file=file.path(files_intermediate, treeGTR.file)) 

