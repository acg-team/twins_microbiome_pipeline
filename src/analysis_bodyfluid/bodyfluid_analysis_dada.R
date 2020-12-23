# @AlexY, created Jan 2020
# analysis of Body Fluid dataset on dada2 ASV
library(adaptiveGPCA)
library(dplyr)
library("RColorBrewer")
library(randomcoloR)
theme_set((theme_bw()))
library(plotly)

source("src/load.R")


project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)

# Load the resulting pipeline file
files_intermediate_dada <- "~/Projects_R/twins_microbiome_pipeline/data_set_bodyfl/files_intermediate_dada"
calculated_ps_file <- "run_BFL_DADA2_Q2_mEE24_trL00_trR00_truncLn210_220_msa_DECIPHER_tax_mapseq.RData"
load(file=file.path(files_intermediate_dada, calculated_ps_file))


#### init: load packages and set path
source("src/configure.R")




############ LOAD DATA and SANITY CHECK
load(file=file.path(metadata_path, "metadata.RData"))
head(df.metadata) # check metadata

# rename the phyloseq object
ps.bfluid <- ps
ps.bfluid

filter.log <- filter.log[-1,]  # delete the first NA row (added by mistake)
colMeans(filter.log)

# transpose OTU matrix if neesesary to make it [taxa x samples] - not neccesary!
if(!taxa_are_rows(ps.bfluid)){
  otu_table(ps.bfluid) <- t(otu_table(ps.bfluid))
  otu_table(ps.bfluid) <- otu_table(ps.bfluid, taxa_are_rows = TRUE)
  print('OTU has been transposed')
}


#### SANITY checks
# check if every sample has at least one taxa and no NA
# TODO: substitute NA with " Not Detected"
sample_sums(ps.bfluid)
get_taxa_unique(ps.bfluid, "Genus")


# OTU: visually check that we have a community matrix: rows are samples (46) and colums are taxa (spesies)
OTU = as(otu_table(ps.bfluid), "matrix")
colnames(OTU) <- c()
dim(OTU)
#View(OTU)

# check the distribution of abandamcies 
hist(OTU[1,],breaks=30)

# TREE: check the validity of tree
tree <- phyloseq::phy_tree(ps.bfluid)
ape::checkValidPhylo(tree)
ape::is.rooted(tree)

# check validity of tree, there should be only 2 children!
edges = tree$edge
mycounts = table(edges[,1]) 
length(mycounts[mycounts == 2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts != 2]) # Number of nodes with more or fewer children
mycounts[mycounts != 2]
#https://github.com/joey711/phyloseq/issues/1215
if(length(mycounts[mycounts != 2])>0){
  print('Your tree has several nodes with more then 2 children! Fixing...')
  # TODO: not sure if that is right approach.. it produces a tree with internal nodes != tips-2
  phy_tree(ps.bfluid) <- ape::multi2di(tree)
}

#plot(phy_tree(tree), show.node.label = TRUE)


# TODO - Interactive choice of r  - https://cran.r-project.org/web/packages/adaptiveGPCA/vignettes/adaptive_gpca_vignette.html


### use all_genera , top30 or all except of top 30 i the analysis
# ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# calculate the total abandance of each genera across all population
genera.sum = tapply(
  X = taxa_sums(ps.bfluid), # sum up all abandancies of all taxa
  INDEX=tax_table(ps.bfluid)[, "Genus"],  # use Genus as an index to sum: 111122223333
  FUN=sum, 
  na.rm=TRUE
  )

n <- length(genera.sum)

all.genera.names <- names(sort(genera.sum, TRUE))
top.genera.names <- all.genera.names[0:30]
no_top.genera.names <- all.genera.names[10:n]


################ START ANALYSIS
# use either all data or for top/ bottom 30
genera.names.to.use <- no_top.genera.names

# TODO: why if I prune for all names I still have less variants? 3209 intead of 4381
physeq.top <- prune_taxa((tax_table(ps.bfluid)[, "Genus"] %in% genera.names.to.use), ps.bfluid)

#physeq.top <- ps.bfluid


# check: the number of taxa should be the same   
ps.bfluid
physeq.top
get_taxa_unique(ps.bfluid, "Genus")

# sanity check: OTU dimensions
otu.top <- as(otu_table(physeq.top), "matrix")
colnames(otu.top) <- c()  # clean the long column names
dim(otu.top)
#View(otu.top)

# each sample has to have non zero variants
sample_sums(physeq.top)
if(any(sample_sums(physeq.top)==0)){
  stop(" one of the samples has zero abandancy for all variants, please check ")
}  

#get_taxa_unique(physeq.top, "Genus")
sample_data(physeq.top)
tax_table(physeq.top)
phy_tree(physeq.top)

########### normalize and log abundancies of community matrix
# do normalization and log transform
# do rel first, only then log - why?!

physeq.top.log <- phyloseq::transform_sample_counts(physeq.top, function(x) log(1 + x)) # log transform with pseudocount
sample_sums(physeq.top.log)

physeq.top.rel <- phyloseq::transform_sample_counts(physeq.top, function(x) x / sum(x) ) # Total Sum Scaling (TSS)
sample_sums(physeq.top.rel)



physeq.top.log.rel <- phyloseq::transform_sample_counts(physeq.top.log, function(x) x / sum(x) ) # Total Sum Scaling (TSS)
sample_sums(physeq.top.log.rel)

physeq.top.rel.log <- phyloseq::transform_sample_counts(physeq.top.rel, function(x) log(1 + x)) # log transform with pseudocount
sample_sums(physeq.top.log)


######
physeq <- physeq.top.rel 


###################  PCoA with weigthed unifrac
# https://joey711.github.io/phyloseq/plot_ordination-examples
dst <- phyloseq::distance(physeq, method="wunifrac", type="samples")
sum(dst==0)   # check if we have NA distances 
# possible warning: https://github.com/joey711/phyloseq/issues/936


# calculate and plot "PCoA" with "unifrac" distances
bf.ord <- phyloseq::ordinate(physeq, "PCoA", "unifrac", weighted=TRUE)
p1 <- phyloseq::plot_ordination(physeq, bf.ord, type="samples", color='Body_site') + geom_density_2d()
print(p1)

# plot in 3D : bf.ord$vectors[,1:3]
plot_ly(x=bf.ord$vectors[,1], y=bf.ord$vectors[,2], z=bf.ord$vectors[,3], 
        type="scatter3d", mode="markers", 
        color=sample_data(physeq)$Body_site, symbol=sample_data(physeq)$State)

# TODO: study outliers?





################## Adaptive gPCA - ordination plot to see distances BTW ALL SAMPLES!
# good source too - https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/doc/MicrobiomeWorkshopII.html
# http://jfukuyama.github.io/adaptiveGPCA/
# extract nessesary matrices for gPCA
pp = adaptiveGPCA::processPhyloseq(physeq)
any(is.na(pp$X))  # if any NA, gPCA wont work
any(is.na(pp$Q))

# run agPCA
out.agpca = adaptiveGPCA::adaptivegpca(pp$X, pp$Q, k = 2)
out.agpca

# PLOT variance explained
plot(out.agpca) 


#### PLOT with color schemas
# Body site colours
bs.colours <- c("red","blue","cyan","burlywood1","magenta")
names (bs.colours) <- levels(df.metadata$Body_site)
names(bs.colours)

# Colour-shape test for phylum - genus
df.tax <- as.data.frame(tax_table(physeq))
df.tax.reduced <- df.tax[,c("Phylum","Genus")]

# package: ‘dplyr’ is a fairly new (2014) package that tries to provide easy tools for the most common data manipulation tasks.
# This addresses a common problem with R in that all operations are conducted in memory and thus the amount of data you can work with is limited by available memory.
# why we use it here? is our dataset big? it is primarily for big data...
df.tax.reduced %>%
  group_by(Phylum) %>%
  summarize(n_unique = n_distinct(Genus))



# TODO: plot exposed / non exposed separatelly
BS <- ggplot(data.frame(out.agpca$U, sample_data(physeq))) +
      geom_point(aes(x = Axis1, y = Axis2, color = Body_site, shape = State)) +
      #geom_density_2d(aes(x = Axis1, y = Axis2, color = Body_site, shape = State)) +  
      ggtitle(paste("Adapt gPCA, sample data ", calculated_ps_file)) +
      #geom_text( aes(label=sample_data(physeq)$SampleID), hjust=0, vjust=0) +
      scale_color_manual(values=bs.colours) +
      xlab("Axis 1") + ylab("Axis 2")
plot(BS)




########## PLOT by taxes
# Note that, it’s possible to plot variables and to color them according to either 
# i) their quality on the factor map (cos2) or 
# ii) their contribution values to the principal components (contrib). 
# colour = Genus

phyla.colours <- brewer.pal(n = 8, name = "Dark2")
names(phyla.colours) <- unique(df.tax.reduced$Phylum)

palette <- distinctColorPalette(n)
#pie(rep(1, n), col=palette)

PHYL <- ggplot(data.frame(out.agpca$QV, tax_table(physeq))) +
  ggtitle(paste("Adaptive gPCA for ", calculated_ps_file)) +
  geom_point(aes(x = Axis1, y = Axis2, shape=Phylum, colour = Genus)) +
  xlab("Pr Axis 1") + ylab("Pr Axis 2") + 
  scale_color_manual(values=palette) #+
  #theme (legend.position ="none")
plot(PHYL)

#t = inspectTaxonomy(out.agpca, physeq)








