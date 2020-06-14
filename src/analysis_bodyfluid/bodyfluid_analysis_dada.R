# @AlexY, created Jan 2020
# analysis of Body Fluid dataset on dada2 ASV


#################### 0: CONFIGURATION : PLEASE SET!  ##############
conf <- vector(mode="list", length=3)
names(conf) <- c("location", "dataset", "pipeline")

### now we set it as calculare on local macbook and use 34 new dataset
conf$location <- "LOCAL"  # LOCAL / HOMESERVER  / ETHSERVER
conf$dataset <- "BODYFL"    #   TWIN / "BODYFL" /
conf$pipeline <- "DADA2"   # QIIME / DADA2

################# FILTERING parameters 
dada_param <- vector(mode="list", length=2)
names(dada_param) <- c("QUALITY_THRESHOLD", "maxEE")
dada_param$QUALITY_THRESHOLD <- 2
dada_param$maxEE <- c(2,4)
dada_param$MSA_aligner <- "MUSCLE"   # DECIPHER  MUSCLE  clustalw 
dada_param$tree_method <- "RAXML"    # PHANGORN   


#### init: load packages and set path
# NOTE - swith from TWIN to BDFLUID!!! in configure.R
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)

source("src/load.R")
source("src/configure.R")

library(adaptiveGPCA)
library(dplyr)
library("RColorBrewer")
library(randomcoloR)
theme_set((theme_bw()))



############ LOAD Budy Fluid PhyloSeq file + metadata
calculated_ps_file <- "phyloseq_BODYFL_DADA2_Q2_maxEE24.RData"
###############

load(file=file.path(metadata_path, "metadata.RData"))
load(file=file.path(files_intermediate_dada, calculated_ps_file))
head(df.metadata) # check metadata

# rename the phyloseq object
ps.bfluid <- ps
ps.bfluid

# transpose OTU matrix if neesesary to make it [taxa x samples] - not neccesary!
if(!taxa_are_rows(ps.bfluid)){
  otu_table(ps.bfluid) <- t(otu_table(ps.bfluid))
  otu_table(ps.bfluid) <- otu_table(ps.bfluid, taxa_are_rows = TRUE)
}



############ SANITY
# check if every sample has at least one taxa and no NA
sample_sums(ps.bfluid)
get_taxa_unique(ps.bfluid, "Phylum") #("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
get_taxa_unique(ps.bfluid, "Class")


# visually check that we have a community matrix: rows are samples (46) and colums are taxa (spesies)
OTU = as(otu_table(ps.bfluid), "matrix")
colnames(OTU) <- c()
dim(OTU)
#View(OTU)

# check the dictribution of abandamcies 
hist(OTU[1,],breaks=20)

# check the tree, thi sis ape package
tree <- phy_tree(ps.bfluid)
ape::checkValidPhylo(tree)
is.rooted(tree)
#plot(phy_tree(tree), show.node.label = TRUE)


# TODO - explore all sequnces names and see NNN, length etc artefacts
# TODO - plot ggplot in 3D axes
# TODO - Interactive choice of r  - https://cran.r-project.org/web/packages/adaptiveGPCA/vignettes/adaptive_gpca_vignette.html


######### Get top 30 genera (need to do it before normalization!) otherwise taxa_sum() does not work
# ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
n <- 30

# total number of each Genera
genera.sum = tapply(taxa_sums(ps.bfluid), tax_table(ps.bfluid)[, "Genus"], sum, na.rm=TRUE)
topGenera = names(sort(genera.sum, TRUE))[1:n]  # names of top30 genera

physeq.top <- prune_taxa((tax_table(ps.bfluid)[, "Genus"] %in% topGenera), ps.bfluid)
otu.top <- as(otu_table(physeq.top), "matrix")
colnames(otu.top) <- c()  # clean the long column names
dim(otu.top)
#View(otu.top)

sample_sums(physeq.top) # sanity again
#get_taxa_unique(physeq.top, "Genus")

# remove taxa with abandancy less then 10
physeq.top = prune_samples(sample_sums(physeq.top)>=10, physeq.top)
sample_sums(physeq.top)


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

###########  Ordination with phyloseq
# https://joey711.github.io/phyloseq/plot_ordination-examples
dst <- phyloseq::distance(physeq, method="wunifrac", type="samples")
sum(dst==0)   # check if we have NA distances 

# calculate and plot "PCoA" with "unifrac" distances
bf.ord <- phyloseq::ordinate(physeq, "PCoA", "unifrac", weighted=TRUE)
p1 <- phyloseq::plot_ordination(physeq, bf.ord, type="samples", color='Body_site')
print(p1)



############# Adaptive gPCA - ordination plot to see distances BTW ALL SAMPLES!
# good source too - https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/doc/MicrobiomeWorkshopII.html
# http://jfukuyama.github.io/adaptiveGPCA/
# extract nessesary matrices for gPCA
pp = adaptiveGPCA::processPhyloseq(physeq)
any(is.na(pp$X))  # if any NA, gPCA wont work
any(is.na(pp$Q))

# run agPCA
out.agpca = adaptiveGPCA::adaptivegpca(pp$X, pp$Q, k = 2)
out.agpca
#tax30 <- inspectTaxonomy(out.agpca, physeq30.rel)


### PLOT results (black and white)
#plot(out.agpca) # variance explained
#plot(out.agpca, type = "samples", axes = c(1,2))   # all samples in 2D, no colod
#plot(out.agpca, type = "variables", axes = c(1,2))  # ??


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


phyla.colours <- brewer.pal(n = 8, name = "Dark2")
names(phyla.colours) <- unique(df.tax.reduced$Phylum)

BS <- ggplot(data.frame(out.agpca$U, sample_data(physeq))) +
      ggtitle(paste("Adap gPCA, sample data ", calculated_ps_file)) +
      geom_point(aes(x = Axis1, y = Axis2, color = Body_site, shape = State)) +
      #geom_text( aes(label=sample_data(physeq30.rel)$SampleID), hjust=0, vjust=0) +
      scale_color_manual(values=bs.colours) +
      xlab("Axis 1") + ylab("Axis 2")
plot(BS)




########## PLOT by taxes
# Note that, it’s possible to plot variables and to color them according to either 
# i) their quality on the factor map (cos2) or 
# ii) their contribution values to the principal components (contrib). 
# colour = Genus

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








