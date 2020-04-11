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
files_intermediate  <- file.path(project_path, "data_set_bodyfl/files_intermediate")
source("src/load.R")
source("src/configure.R")
theme_set((theme_bw()))


############ LOAD Budy Fluid PhyloSeq file + metadata
metadata_path   <- file.path(project_path, "data_set_bodyfl/metadata")
load(file=file.path(metadata_path, "metadata.RData"))
load(file=file.path(files_intermediate, "phyloseq_object.RData"))
head(df.metadata)

# rename the phyloseq object
ps.bfluid <- ps.tweens

# sanity check of loaded data
get_taxa_unique(ps.bfluid, "Phylum")
get_taxa_unique(ps.bfluid, "Genus")
# Q: is NA a problem?
sample_data(ps.bfluid)

# visually check that we hava a community matrix: rows are samples (46) and colums are taxa (spesies)
a <- ps.bfluid@otu_table@.Data
colnames(a) <- c()


########### normalize and log abundancies of community matrix
# do normalization and log transform
ps.bfluid.log <- phyloseq::transform_sample_counts(ps.bfluid, function(x) log(1 + x)) # log transform with pseudocount
any(is.na(ps.bfluid.log@otu_table@.Data))

ps.bfluid.rel <- phyloseq::transform_sample_counts(ps.bfluid, function(x) x / sum(x) ) # Total Sum Scaling (TSS)
any(is.na(ps.bfluid.rel@otu_table@.Data))  # have NAs here! cannot use!


# Which ps to use?
ps <- ps.bfluid.log

# WHY ??
# ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
# Top 30 genera - TSS without blood
genera.sum = tapply(taxa_sums(ps), tax_table(ps)[, "Genus"], sum, na.rm=TRUE)
top30Genera = names(sort(genera.sum, TRUE))[1:30]
physeq30 <- prune_taxa((tax_table(ps)[, "Genus"] %in% top30Genera), ps)

# check that we haev only 30 genus
get_taxa_unique(physeq30, "Genus")



########## do Adaptive gPCA - ordination plot to see distances BTW ALL SAMPLES!
# on LOG abandancies

# extract nessesary matrices for gPCA
pp = adaptiveGPCA::processPhyloseq(physeq30)
# if any NA, gPCA wont work
any(is.na(pp$X))
any(is.na(pp$Q))

# run agPCA
out.agpca = adaptiveGPCA::adaptivegpca(pp$X, pp$Q, k = 2)
out.agpca
#tax30 <- inspectTaxonomy(out.agpca, ps.bfluid)


################## PLOT 
plot(out.agpca) # variance explained
plot(out.agpca, type = "samples", axes = c(1,2))   # all samples in 2D, no colod
plot(out.agpca, type = "variables", axes = c(1,2))  # ??

########### PLOT with color schemas
# Body site colours
bs.colours <- c("red","blue","cyan","burlywood1","magenta")
#metadata_nb <- df.metadata[!df.metadata$Body_site=="blood"]   # blood is already deleted
names (bs.colours) <- levels(df.metadata$Body_site)
names(bs.colours)

# Colour-shape test for phylum - genus
df.tax <- as.data.frame(tax_table(physeq30))
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
  data.frame(out.agpca$U, sample_data(physeq30))) +
  geom_point(aes(x = Axis1, y = Axis2, color = Body_site, shape = State)) +
  scale_color_manual(values=bs.colours) +
  xlab("Axis 1") + ylab("Axis 2")
plot(BS)

# compare it with results of QIIME?



########## WHAT is THAT?
# colour = Genus
n <- 30
palette <- distinctColorPalette(n)
pie(rep(1, n), col=palette)

PHYL <- ggplot(
  data.frame(out.agpca$QV, tax_table(physeq30))) +
  geom_point(aes(x = Axis1, y = Axis2, shape=Phylum, colour = Genus)) +
  xlab("Axis 1") + ylab("Axis 2") + 
  scale_color_manual(values=palette)+
  theme (legend.position ="none")
plot(PHYL)


############ alpha diversity

