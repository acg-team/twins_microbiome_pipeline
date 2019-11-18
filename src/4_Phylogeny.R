#### INFER A PHYLOGENY TREE
##############################################
# http://www.metagenomics.wiki/tools/phylogenetic-tree


#### init: load packages and set path
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)
source("src/load_initialize.R")

load(file=file.path(files_intermediate, seqtab.file)) 
load(file=file.path(files_intermediate, seqtab.snames.file)) 
load(file=file.path(files_intermediate, taxtab.file))


##############  MSA Construction
seqs <- dada2::getSequences(seqtab)   # does not work for some reason
names(seqs) <- seqtab.samples.names # This propagates to the tip labels of the tree

#msa package provides a unified R/Bioconductor interface to MSA (ClustalW, ClustalOmega, Muscle)
# TODO: check for ClustalW specific parameters
#microbiome.msa.clustalW <- msa::msa(seqs, method="ClustalW", type="dna", order="input")

# TODO: check for Muscle specific parameters
#microbiome.msa.muscle <- msa::msa(seqs, method="Muscle", type="dna", order="input")
tic()
microbiome.msa.muscle <- msa::msaMuscle(seqs, type="dna", order="input")

# TODO:  visualize MSA , type: msa::MsaDNAMultipleAlignment
msa::msaPrettyPrint(x=microbiome.msa.muscle, output="pdf", subset=NULL,file=paste0("msa.muscle", ".pdf"))
               
# save MSA to a file 
save(microbiome.msa.clustalW,microbiome.msa.muscle, file=file.path(files_intermediate, msa.file)) 




########## Construct tree 

################# NJ tree with phangorn
tic()
# infer a guide tree
phang.align <- phangorn::as.phyDat(mult, type="DNA", names=seqtab.samples.names)
dm <- dist.ml(phang.align)
treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
toc()

################ ML tree with phangorn
## WARNING! Might take a lot of time
# try ML tree with Jukes-Cantor model
tic()
fit = phangorn::pml(treeNJ, data=phang.align)
toc()

# try ML tree with GTR model
tic()
fitGTR <- update(fit, k=4, inv=0.2)  #???
fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
toc()

############### ML tree with RAxML: ML tree for species >1000 with fast heuristics
# NOTE: need to install raxml on local MAC first
exec.path <- "/Users/alex/bioinf_tools/RAxML/raxmlHPC-PTHREADS-AVX"
exec.path.ubuntu <- "/home/alex/installed/BIOINF_tools/RAxML/raxmlHPC-PTHREADS-AVX"

# msa data must be in DNAbin format??
msa.raxm <- microbiome.msa.muscle

tic()
# f - RAxML algorithm
# N - Integers give the number of independent searches on different starting tree or replicates in bootstrapping. 
# p - Integer, setting a random seed for the parsimony starting trees.
# return tr is a list of tr[1] - info, tr[2] - best tree 
tr <- raxml(msa.raxm, m = "GTRGAMMA", f = "d", N = 1, p = 1234, exec = exec.path.ubuntu, threads=2, file="twin_tree") 
toc()


#########################

# TODO: save trees in appropriate format

# save the tree to file
save(fitGTR, file=file.path(files_intermediate, treeGTR.file)) 

