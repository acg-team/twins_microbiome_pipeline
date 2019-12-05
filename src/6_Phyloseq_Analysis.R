### use Phylogeq oblest (ps.tweens) 

#### init: load packages and set path
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)
source("src/load_initialize.R")

### LOAD PREVIOUS DATA
load(file=file.path(files_intermediate, phyloseq.file)) 


#########  Start exploration and analysis
# https://joey711.github.io/phyloseq-demo/unifrac.html
# do we need a rooted tree?

# Unifrac assess a distance btw two sets of microbial community based on their tree position and abandunce 
# this might take some time!
# 7 hours
tic()
unifrac.dist.matr <- phyloseq::distance(ps.tweens, method="unifrac", type="samples", fast=TRUE, parallel=TRUE)
toc()
# NOTE: Randomly assigning root as 

save(unifrac.dist.matr, file=file.path(files_intermediate, phyloseq_analysis.file)) 





# TODO: fither analysis of sample by sample UNIFRAC distance matrix?
# clustering? PCoA? PCA? 
# PCoA is recommended over PCA when there are lots of missing data and when there are fewer individuals than characters




# TODO figure out what genera is most usefull
# filter common taxa
# check if microbiota is stable thouth the life




##################### TODO: not finished 
# melt all sequences of one taxa to one abundance
# https://github.com/joey711/phyloseq/issues/418
glom <- phyloseq::tax_glom(ps, taxrank = 'Family')   # aglomerate taxa
dat <- psmelt(glom)   # convert to df
dat$Family <- as.character(dat$Family)  # convert Phylum to a character vector from a factor

Phylum_abundance <- aggregate(Abundance~Sample+Phylum, dat, FUN=sum)  #aggregate

# reorganize the table so that each phylum is a column
Phylum_abundance <- cast(Phylum_abundance, Sample ~ Phylum)

save(Phylum_abundance, file=file.path(result_path, "Phylum_abundance.RData")) 

write.table(Phylum_abundance, file.path(result_path,"Famuly_abundance.csv"), sep=",")