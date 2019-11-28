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
seqs <- dada2::getSequences(seqtab)   # 8299
names(seqs) <- seqs    # This propagates to the tip labels of the tree
seq.variant.names <- names(seqs)
seq.number <- length(seq.variant.names)
# generate short names (RAXML requires names to be less then 256)
prefix <- "sv_seq_variant"
suffix <- seq(1:seq.number)
seq.variant.short.names <- paste(prefix, suffix, sep='_')
names(seqs) <- seq.variant.short.names


#msa package provides a unified R/Bioconductor interface to MSA (ClustalW, ClustalOmega, Muscle)
# TODO: check for ClustalW specific parameters
#microbiome.msa.clustalW <- msa::msa(seqs, method="ClustalW", type="dna", order="input")

# TODO: look for another package for MSA, this on msa is very badly written
# TODO: check for Muscle specific parameters
tic()
microbiome.msa.muscle <- msa::msaMuscle(seqs, type="dna", order="input")
print("msa (muscle) took:")
toc()  # 5 hours
print(microbiome.msa.muscle)
writeXStringSet(unmasked(microbiome.msa.muscle), file=file.path(result_path, "msa_muscle.fasta"))

tic()
microbiome.msa.clustalw <- msa::msaClustalW(seqs, type="dna", order="input")
print("msa (clustalw) took:")
toc() # 6 hours
print(microbiome.msa.clustalw)
#microbiome.msa.clustalw@unmasked@ranges@NAMES[3000:4000]
writeXStringSet(unmasked(microbiome.msa.clustalw), file=file.path(result_path, "msa_clustalw.fasta"))

# save MSA to a file 
save(microbiome.msa.muscle, microbiome.msa.clustalw, file=file.path(files_intermediate, msa.file)) 

# TODO:  visualize MSA , type: msa::MsaDNAMultipleAlignment or use UGene browser
# msa::msaPrettyPrint(x=microbiome.msa.muscle, output="tex", subset=NULL)
# tools::texi2pdf("msaPrettyPrintOutput.tex",clean=TRUE)
              



########## Construct tree 

################# NJ tree with phangorn
my.msa <- microbiome.msa.muscle

tic()
# infer a guide tree
phang.align <- phangorn::as.phyDat(my.msa, type="DNA", names=seqtab.samples.names)
dm <- dist.ml(phang.align)
treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
toc()
save(treeNJ, file=file.path(files_intermediate, phylo.file)) 

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

# save the tree to file
save(treeNJ, fitGTR, file=file.path(files_intermediate, phylo.file)) 



############### ML tree with RAxML: ML tree for species >1000 with fast heuristics
# NOTE: need to install raxml on local MAC first
exec.path <- "/Users/alex/bioinf_tools/RAxML/raxmlHPC-PTHREADS-AVX"
exec.path.ubuntu <- "/home/alex/installed/BIOINF_tools/RAxML/raxmlHPC-PTHREADS-AVX"

# convert msa::MsaDNAMultipleAlignment data into ips::DNAbin (ape::DNAbim) format!

msa.dnabin <- msa::msaConvert(my.msa, "ape::DNAbin")

# vizual control of MSA
labels(msa.dnabin)
print(msa.dnabin)

save(msa.dnabin, file=file.path(files_intermediate, msa.file)) 

# I exported my alignment to a server to run RAXML, it took 9 days for 5000 sequences.
# f - RAxML algorithm
# N - Integers give the number of independent searches on different starting tree or replicates in bootstrapping. 
# p - Integer, setting a random seed for the parsimony starting trees.
# return tr is a list of tr[1] - info, tr[2] - best tree 
tic()
tr <- raxml(msa.dnabin, f = "d", N = 2, p = 1234, exec = exec.path.ubuntu, threads=4) # , file="RAxMLtwin_tree",  m = "GTRGAMMA",
toc()


#########################

# TODO: save trees in appropriate format



