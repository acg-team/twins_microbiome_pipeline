#### do MSA and INFER A PHYLOGENY TREE with
# MSA - a core Bioconductor package, for multiple sequence alignment
# Phangorn (NJ, ML Felsenstein)
# RAxML (TODO)
# ape - Analyses of Phylogenetics and Evolution
##############################################
# http://www.metagenomics.wiki/tools/phylogenetic-tree


#### init: load packages and set path
load(file=file.path(files_intermediate_dada, seqtab.file)) 
load(file=file.path(files_intermediate_dada, seqtab.snames.file)) 


##############  MSA Construction ##############
# extract DNA seq from seqtab object
seqs <- asv_sequences  # NOTE: names of this vector will propagate to the tip labels of the tree


# note: msa package provides a unified R/Bioconductor interface to MSA (ClustalW, ClustalOmega, Muscle)
# TODO: check for ClustalW specific parameters
#microbiome.msa.clustalW <- msa::msa(seqs, method="ClustalW", type="dna", order="input")

# TODO: look for another package for MSA, this on msa is very badly written
# TODO: check for Muscle specific parameters

# option 0:  AlignSeqsfrom the DECIPHER
if (tools_param$MSA_aligner=="DECIPHER"){
  print("--> run MSA by DECIPHER")
  microbiome.msa.decipher <- DECIPHER::AlignSeqs( DNAStringSet(seqs) )
  Biostrings::writeXStringSet(microbiome.msa.decipher, file=file.path(result_path, "msa_decipher.fasta"))
  save(microbiome.msa.decipher, file=file.path(files_intermediate_dada, msa.file)) 
  }


# option 1: generate MSA with Muscle
# 5 hours
if (tools_param$MSA_aligner=="MUSCLE"){
  print("--> run MSA by MUSCLE")
  tic()
  microbiome.msa.muscle <- msa::msaMuscle(seqs, type="dna", order="input")
  cat("msa (muscle) took: ")
  toc()  # 1972.561sec
  print(microbiome.msa.muscle)
  rownames(microbiome.msa.muscle)
  
  # save MSA as a fasta file for possible vizualization with UGene browser
  Biostrings::writeXStringSet(unmasked(microbiome.msa.muscle), file=file.path(result_path, "msa_muscle.fasta"))
  save(microbiome.msa.muscle, file=file.path(files_intermediate_dada, msa.file)) 
}



# option 2: generate MSA with clustalW
# 6 hours
if (tools_param$MSA_aligner=="clustalw"){
  print("--> run MSA by clustalw")
  tic()
  microbiome.msa.clustalw <- msa::msaClustalW(seqs, type="dna", order="input")
  cat("msa (clustalw) took: ")
  toc() 
  print(microbiome.msa.clustalw)
  Biostrings::writeXStringSet(unmasked(microbiome.msa.clustalw), file=file.path(result_path, "msa_clustalw.fasta"))
  
  # save objects for reusing late in pipeline 
  save(microbiome.msa.clustalw,  file=file.path(files_intermediate_dada, msa.file)) 
}

              

#################################################
# use one of this 
if (tools_param$MSA_aligner=="DECIPHER"){ my.msa <- microbiome.msa.decipher }
if (tools_param$MSA_aligner=="MUSCLE"){ my.msa <- microbiome.msa.muscle }
if (tools_param$MSA_aligner=="clustalw"){ my.msa <- microbiome.msa.clustalw }




####################### Infer a phylogenetic tree 
# TODO: use a separate file for each method (like taxonomy)

########### OPTION 1:  fast NJ tree, can be used as guide tree as well

# infer a tree with fast NJ method 
# 40 min
if (tools_param$tree_method=="PHANGORN"){
  print("-------> Tree inference by PHANGORN")
  tic()
  phang.align <- phangorn::as.phyDat(my.msa, type="DNA", names=samples.names)
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
  fitJC <- phangorn::pml(tree=treeNJ, data=phang.align)   # pmlcomputes  the  likelihood  of  a  phylogenetic  tree 
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

save(my.msa, msa.dnabin, file=file.path(files_intermediate_dada, msa.file)) 

if (tools_param$tree_method=="RAXML"){
  print("-------> Tree inference by RAXML")
  # Parameters:
  # f - RAxML algorithm
  # N - Integers give the number of independent searches on different starting tree or replicates in bootstrapping. 
  # p - Integer, setting a random seed for the parsimony starting trees.
  # return tr is a list of tr[1] - info, tr[2] - best tree (rooted)
  
  # 5.5h
  tic()
  tree.raxml <- ips::raxml(
    as.matrix(msa.dnabin), 
    m = "GTRGAMMA",
    f = "d",   # d - new rapid hill-climbing / "a", # best tree and bootstrap
    N = 4, # number of bootstrap replicates
    p = 2234, # random number seed
    exec = raxm.exec.path, 
    threads=6,
    file="RAxML_tree"
  ) # , file="RAxMLtwin_tree",  m = "GTRGAMMA",
  
  cat("RAXML elapsed time: ")
  toc() # 3045.909 sec = 0.8 h om 6 core server - very fast
  
  my.tree <- tree.raxml
  save(my.tree, file=file.path(files_intermediate_dada, phylo.file)) 
  
}







