# @AlexY, created Jan 2020
# do analysis of Body Fluid dataset

library(adaptiveGPCA)
library(dplyr)
library("RColorBrewer")
library(randomcoloR)
theme_set((theme_bw()))

#### init: load packages and set path
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)
files_intermediate  <- file.path(project_path, "data_set_bodyfl/files_intermediate")
source("src/load.R")
source("src/configure.R")


### LOAD Budy Fluid PhyloSeq file + metadata
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

a<-ps.bfluid@otu_table@.Data
colnames(a) <- c()


# do normalization and log transform
ps.bfluid.log <- phyloseq::transform_sample_counts(ps.bfluid, function(x) log(1 + x)) # log transform with pseudocount
ps.bfluid.rel <- phyloseq::transform_sample_counts(ps.bfluid, function(x) x / sum(x) ) # Total Sum Scaling (TSS)

#Top 30 genera - TSS without blood
genera.sum = tapply(taxa_sums(ps.bfluid.rel), tax_table(ps.bfluid.rel)[, "Genus"], sum, na.rm=TRUE)
top30Genera = names(sort(genera.sum, TRUE))[1:30]
physeq30_rel = prune_taxa((tax_table(ps.bfluid.rel)[, "Genus"] %in% top30Genera), ps.bfluid.rel)
physeq30_rel_logs <-transform_sample_counts(physeq30_rel, function(x) log(0.00001+x))


#### do Adaptive gPCA

# 2 - processPhyloseq, create a loist of nessesary matrix for gPCA
pp = processPhyloseq(physeq30_rel_logs)


# 3 - run and plot
out.agpca = adaptivegpca(pp$X, pp$Q, k = 2)
out.agpca
#tax30 <- inspectTaxonomy(out.agpca, ps.bfluid)

plot(out.agpca) #defaul is scree
plot(out.agpca, type = "samples", axes = c(1,2))
plot(out.agpca, type = "variables", axes = c(1,2))

# 4 - add color schemas
# Body site colours
bs.colours <- c("red","blue","cyan","burlywood1","magenta")
#metadata_nb <- df.metadata[!df.metadata$Body_site=="blood"]   # why?!
names (bs.colours) <- levels(df.metadata$Body_site)
names(bs.colours)

# Colour-shape test for phylum - genus
df.tax <- as.data.frame(tax_table(physeq30_rel_logs))
df.tax.reduced <- df.tax[,c("Phylum","Genus")]

# package: ‘dplyr’ is a fairly new (2014) package that tries to provide easy tools for the most common data manipulation tasks.
# This addresses a common problem with R in that all operations are conducted in memory and thus the amount of data you can work with is limited by available memory.
# why we use it here? is our dataset big? it is primarily for big data...
df.tax.reduced %>%
  group_by(Phylum) %>%
  summarize(n_unique = n_distinct(Genus))


phyla.colours <- brewer.pal(n = 8, name = "Dark2")
names(phyla.colours) <- unique(df.tax.reduced$Phylum)

BS <- ggplot(data.frame(out.agpca$U, sample_data(physeq30_rel_logs))) +
  geom_point(aes(x = Axis1, y = Axis2, color = Body_site, shape = State)) +
  scale_color_manual(values=bs.colours) +
  xlab("Axis 1") + ylab("Axis 2")


# data.frame(out.agpca$QV, tax_table(physeq30_rel_logs)) : arguments imply differing number of rows: 801, 738
n <- 30
palette <- distinctColorPalette(n)
pie(rep(1, n), col=palette)

PHYL <- ggplot(data.frame(out.agpca$QV, tax_table(physeq30_rel_logs))) +
  geom_point(aes(x = Axis1, y = Axis2, shape= Phylum, colour = Genus)) +
  xlab("Axis 1") + ylab("Axis 2") + 
  scale_color_manual(values=palette)+
  theme (legend.position ="none")
plot(PHYL)


