#### do MSA and INFER A PHYLOGENY TREE with
# MSA - a core Bioconductor package, for multiple sequence alignment
# Phangorn (NJ, ML Felsenstein)
# RAxML (TODO)
# ape - Analyses of Phylogenetics and Evolution
##############################################
# http://www.metagenomics.wiki/tools/phylogenetic-tree


#### init: load packages and set path
source("src/load.R")
source("src/configure.R")
setwd(project_path)

load(file=file.path(files_intermediate_dada, seqtab.file)) 
load(file=file.path(files_intermediate_dada, seqtab.snames.file)) 
load(file=file.path(files_intermediate_dada, taxtab.file))


##############  MSA Construction ##############
# extract DNA seq from seqtab object
seqs <- dada2::getSequences(seqtab)   # 8299
names(seqs) <- seqs    # This propagates to the tip labels of the tree

# define variables
seq.variant.names <- names(seqs)
seq.number <- length(seq.variant.names)

# generate short names (RAXML requires names to be less then 256)
# TODO - a questional... may be use the first 10 letters? with number?
prefix <- "sv_seq_variant"
suffix <- seq(1:seq.number)
seq.variant.short.names <- paste(prefix, suffix, sep='_')
names(seqs) <- seq.variant.short.names


# note: msa package provides a unified R/Bioconductor interface to MSA (ClustalW, ClustalOmega, Muscle)
# TODO: check for ClustalW specific parameters
#microbiome.msa.clustalW <- msa::msa(seqs, method="ClustalW", type="dna", order="input")

# TODO: look for another package for MSA, this on msa is very badly written
# TODO: check for Muscle specific parameters

# option 0:  AlignSeqsfrom the DECIPHER
if (dada_param$MSA_aligner=="DECIPHER"){
  print("--> run MSA by DECIPHER")
  microbiome.msa.decipher <- DECIPHER::AlignSeqs( DNAStringSet(seqs) )
  Biostrings::writeXStringSet(microbiome.msa.decipher, file=file.path(result_path, "msa_decipher.fasta"))
  save(microbiome.msa.decipher, seq.variant.names, file=file.path(files_intermediate_dada, msa.file)) 
  }


# option 1: generate MSA with Muscle
# 5 hours
if (dada_param$MSA_aligner=="MUSCLE"){
  print("--> run MSA by MUSCLE")
  tic()
  microbiome.msa.muscle <- msa::msaMuscle(seqs, type="dna", order="input")
  print("msa (muscle) took:")
  toc()  # 1972.561sec
  print(microbiome.msa.muscle)
  rownames(microbiome.msa.muscle)
  
  # save MSA as a fasta file for possible vizualization with UGene browser
  Biostrings::writeXStringSet(unmasked(microbiome.msa.muscle), file=file.path(result_path, "msa_muscle.fasta"))
  save(microbiome.msa.muscle,seq.variant.names, file=file.path(files_intermediate_dada, msa.file)) 
}



# option 2: generate MSA with clustalW
# 6 hours
if (dada_param$MSA_aligner=="clustalw"){
  print("--> run MSA by clustalw")
  tic()
  microbiome.msa.clustalw <- msa::msaClustalW(seqs, type="dna", order="input")
  print("msa (clustalw) took:")
  toc() 
  print(microbiome.msa.clustalw)
  #microbiome.msa.clustalw@unmasked@ranges@NAMES[3000:4000]
  Biostrings::writeXStringSet(unmasked(microbiome.msa.clustalw), file=file.path(result_path, "msa_clustalw.fasta"))
  
  # save objects for reusing late in pipeline 
  save(microbiome.msa.clustalw, seq.variant.names, file=file.path(files_intermediate_dada, msa.file)) 
}

              

#################################################
# use on of this 
if (dada_param$MSA_aligner=="DECIPHER"){ my.msa <- microbiome.msa.decipher }
if (dada_param$MSA_aligner=="MUSCLE"){ my.msa <- microbiome.msa.muscle }
if (dada_param$MSA_aligner=="clustalw"){ my.msa <- microbiome.msa.clustalw }



####################### Infer a phylogenetic tree 

########### OPTION 1:  fast NJ tree, can be used as guide tree as well

# TODO - choose only one methor of tree
# infer a tree with fast NJ method 
# 40 min
if (dada_param$tree_method=="PHANGORN"){
  print("---> Tree inference by PHANGORN")
  tic()
  phang.align <- phangorn::as.phyDat(my.msa, type="DNA", names=seqtab.samples.names)
  dm <- dist.ml(phang.align)  #distance matrix
  treeNJ <- phangorn::NJ(dm) # "phylo" object (a tree)
  toc()
  save(treeNJ, file=file.path(files_intermediate_dada, phylo.file)) 
  
  ####### refine NJ tree with nt substitution model by Felsenstein ML mehod
  ## WARNING! Might take a lot of time
  # infer ML tree with Jukes-Cantor model (JC69, default one), usin NJ as a guide tree
  # fitJC is "pml" object, tree can be extracted as fitJC$tree, also has logLik etc parameters
  # 55 min
  tic()
  fitJC = phangorn::pml(tree=treeNJ, data=phang.align)   # pmlcomputes  the  likelihood  of  a  phylogenetic  tree 
  fitJC <- optim.pml(fitJC)    # optimize edge length etc parameters
  toc()
  save(treeNJ, fitJC, file=file.path(files_intermediate_dada, phylo.file))
  
  # futher refine ML tree with GTR+G+I model
  tic()
  # change parameters of pml: k=Number of intervals of the discrete gamma distribution, inv=Proportion of invariable sites
  # What is that parameters?!
  fitGTR <- update(fitJC, k=4, inv=0.2)  
  fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
  toc()
  
  # save the tree to file
  save(treeNJ, fitJC, fitGTR, file=file.path(files_intermediate_dada, phylo.file)) 
  
  my.tree <- fitGTR
}



############### OPTION2 :  ML tree with RAxML: ML tree for species >1000 with fast heuristics
####  NOTE: need to install raxml on local MAC first
#exec.path.mac <- "/Users/alex/bioinf_tools/RAxML/raxmlHPC-PTHREADS-AVX"
#exec.path.ubuntu <- "/home/alex/installed/BIOINF_tools/RAxML/raxmlHPC-PTHREADS-AVX"

# convert msa::MsaDNAMultipleAlignment data into ips::DNAbin (ape::DNAbim) format!
msa.dnabin <- ape::as.DNAbin(my.msa, check.names=TRUE)


# vizual control of MSA
labels(msa.dnabin)
print(msa.dnabin)   # Base composition: acgt = NaN! why?

save(my.msa, seq.variant.names, msa.dnabin, file=file.path(files_intermediate_dada, msa.file)) 

if (dada_param$tree_method=="RAXML"){
  # Parameters:
  # f - RAxML algorithm
  # N - Integers give the number of independent searches on different starting tree or replicates in bootstrapping. 
  # p - Integer, setting a random seed for the parsimony starting trees.
  # return tr is a list of tr[1] - info, tr[2] - best tree (rooted)
  

#  alignment.rax.gtr <- raxml(alignment,
#                             m="GTRGAMMAIX", # model
#                             f="a", # best tree and bootstrap
#                             p=1234, # random number seed
#                             x=2345, # random seed for rapid bootstrapping
#                             N=100, # number of bootstrap replicates
#                             file="alignment", # name of output files
#                             exec="raxmlHPC-PTHREADS-SSE3", # name of executable
#                             threads=20
#  )

  
  # 5.5h
  tic()
  tree.raxml <- ips::raxml(
    as.matrix(msa.dnabin), 
    m = "GTRGAMMA",
    f = "d",   # d - new rapid hill-climbing / "a", # best tree and bootstrap
    N = 3, 
    p = 1234, 
    exec = raxm.exec.path, 
    threads=6,
    file="RAxML_tree"
  ) # , file="RAxMLtwin_tree",  m = "GTRGAMMA",
  
  toc() # 3045.909 sec = 0.8 h om 6 core server - very fast
  
  my.tree <- tree.raxml
  save( my.tree, file=file.path(files_intermediate_dada, phylo.file)) 
  
}







