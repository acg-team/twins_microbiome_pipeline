# @AlexY, created Jan 2020
# Alpha and Beta diversity analysis of Body Fluid dataset on dada2 ASV

library(adaptiveGPCA)
library(dplyr)
library("RColorBrewer")
library(randomcoloR)

#### init: load packages and set path
# NOTE - swith from TWIN to BDFLUID!!! in configure.R
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)

source("src/load.R")
source("src/configure.R")
theme_set((theme_bw()))


############ LOAD Budy Fluid PhyloSeq file + metadata
load(file=file.path(metadata_path, "metadata.RData"))
load(file=file.path(files_intermediate_dada, "phyloseq_object.RData"))
head(df.metadata) # check metadata

# rename the phyloseq object
ps.bfluid <- ps.tweens
ps.bfluid

# transpose OTU matrix if neesesary to make it [taxa x samples] - not neccesary!
if(!taxa_are_rows(ps.bfluid)){
  otu_table(ps.bfluid) <- t(otu_table(ps.bfluid))
  otu_table(ps.bfluid) <- otu_table(ps.bfluid, taxa_are_rows = TRUE)
}



#### SANITY
# check if every sample has at least one taxa and no NA
sample_sums(ps.bfluid)
get_taxa_unique(ps.bfluid, "Phylum") #("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
is.rooted(phy_tree(ps.bfluid))

# visually check that we hava a community matrix: rows are samples (46) and colums are taxa (spesies)
OTU = as(otu_table(ps.bfluid), "matrix")
colnames(OTU) <- c()



##### Get top 30 genera (need to do it before normalization!) otherwise taxa_sum() does not work
# ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
n <- 30

genera.sum = tapply(taxa_sums(ps.bfluid), tax_table(ps.bfluid)[, "Genus"], sum, na.rm=TRUE)
top30Genera = names(sort(genera.sum, TRUE))[1:n]
physeq30 <- prune_taxa((tax_table(ps.bfluid)[, "Genus"] %in% top30Genera), ps.bfluid)
otu.30 <- as(otu_table(physeq30), "matrix")

sample_sums(physeq30) # sanity again




########### normalize and log abundancies of community matrix
# do normalization and log transform
# do rel first, only then log

physeq30.rel <- phyloseq::transform_sample_counts(physeq30, function(x) x / sum(x) ) # Total Sum Scaling (TSS)
sample_sums(physeq30.rel)

#physeq30.log <- phyloseq::transform_sample_counts(physeq30, function(x) log(1 + x)) # log transform with pseudocount




###########  Ordination
# NOTE: something wrong with distance matrix
dst <- phyloseq::distance(physeq30.rel, method="wunifrac", type="samples")
sum(dst==0)   # check if we have NA distances 

# calculate and plot "PCoA" with "unifrac" distances
bf.ord <- ordinate(physeq30.rel, "PCoA", "unifrac", weighted=TRUE)
p1 <- plot_ordination(physeq30.rel, bf.ord, type="samples", color='Body_site')
print(p1)



############# Adaptive gPCA - ordination plot to see distances BTW ALL SAMPLES!
# http://jfukuyama.github.io/adaptiveGPCA/
# extract nessesary matrices for gPCA
pp = adaptiveGPCA::processPhyloseq(physeq30.rel)
any(is.na(pp$X))  # if any NA, gPCA wont work
any(is.na(pp$Q))

# run agPCA
out.agpca = adaptiveGPCA::adaptivegpca(pp$X, pp$Q, k = 2)
out.agpca
#tax30 <- inspectTaxonomy(out.agpca, physeq30.rel)


### PLOT results
plot(out.agpca) # variance explained
plot(out.agpca, type = "samples", axes = c(1,2))   # all samples in 2D, no colod
plot(out.agpca, type = "variables", axes = c(1,2))  # ??


#### PLOT with color schemas
# Body site colours
bs.colours <- c("red","blue","cyan","burlywood1","magenta")
names (bs.colours) <- levels(df.metadata$Body_site)
names(bs.colours)

# Colour-shape test for phylum - genus
df.tax <- as.data.frame(tax_table(physeq30.rel))
df.tax.reduced <- df.tax[,c("Phylum","Genus")]

# package: ‘dplyr’ is a fairly new (2014) package that tries to provide easy tools for the most common data manipulation tasks.
# This addresses a common problem with R in that all operations are conducted in memory and thus the amount of data you can work with is limited by available memory.
# why we use it here? is our dataset big? it is primarily for big data...
df.tax.reduced %>%
  group_by(Phylum) %>%
  summarize(n_unique = n_distinct(Genus))


phyla.colours <- brewer.pal(n = 8, name = "Dark2")
names(phyla.colours) <- unique(df.tax.reduced$Phylum)

BS <- ggplot(
  data.frame(out.agpca$U, sample_data(physeq30.rel))) +
  geom_point(aes(x = Axis1, y = Axis2, color = Body_site, shape = State)) +
  #geom_text( aes(label=sample_data(physeq30.rel)$SampleID), hjust=0, vjust=0) +
  scale_color_manual(values=bs.colours) +
  xlab("Axis 1") + ylab("Axis 2")
plot(BS)




########## PLOT by taxes?
# Note that, it’s possible to plot variables and to color them according to either 
# i) their quality on the factor map (cos2) or 
# ii) their contribution values to the principal components (contrib). 
# colour = Genus

palette <- distinctColorPalette(n)
#pie(rep(1, n), col=palette)

PHYL <- ggplot(
  data.frame(out.agpca$QV, tax_table(physeq30.rel))) +
  geom_point(aes(x = Axis1, y = Axis2, shape=Phylum, colour = Genus)) +
  xlab("Axis 1") + ylab("Axis 2") + 
  scale_color_manual(values=palette) #+
  #theme (legend.position ="none")
plot(PHYL)

t = inspectTaxonomy(out.agpca, physeq30.rel)








