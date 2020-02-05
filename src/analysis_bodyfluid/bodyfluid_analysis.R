# @AlexY, created Jan 2020
# do analysis of Body Fluid dataset

library(adaptiveGPCA)
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

#### do Adaptive gPCA

# 1 - do normalization and log transform


# 2 - processPhyloseq, create a loist of nessesary matrix for gPCA
pp = adaptiveGPCA::processPhyloseq(ps.bfluid)

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
metadata_nb <- df.metadata[!df.metadata$Body_site=="blood"]   # why?!
names (bs.colours) <- levels(metadata_nb$Body_site)
names(bs.colours)

library("RColorBrewer")
phyla.colours <- brewer.pal(n = 8, name = "Dark2")
names(phyla.colours) <- unique(df.tax.reduced$Phylum)

BS <- ggplot(data.frame(out.agpca$U, sample_data(physeq30_rel_logs))) +
  geom_point(aes(x = Axis1, y = Axis2, color = Body_site, shape = State)) +
  scale_color_manual(values=bs.colours) +
  xlab("Axis 1") + ylab("Axis 2")
