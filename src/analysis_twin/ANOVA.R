# This is a fule for AKOS to modify
# https://github.com/joey711/phyloseq/wiki/ordinate


source("src/configure.R")
source("src/load.R")
my_path  <- file.path(project_path, "src/analysis_twin")

load(file=file.path(metadata_path, metadata.file))
load(file=file.path(files_intermediate, phyloseq.file))
load('data_set_twin/analysis/unifrac.RData')  # load unifrac distance tables

print(conf)  # check current dataset

# CHECK: here should be a metadata 
df.metadata.4timepoints

# CHECK: here should be a phyloseq object
ps.tweens

# abanduncy normalization and log transform
ps.tweens.norm <- transform_sample_counts(ps.tweens, function(x) x / sum(x) )
ps.tweens.log <- transform_sample_counts(ps.tweens.norm, function(x) log(1 + x))



###############  study within family sample variability   ########
twin.families <- unique(df.metadata.4timepoints$family_id)

print_family <- function(family_id){
  print(df.metadata.4timepoints[df.metadata.4timepoints$family_id==twin.families[family.number], ])
}


# LOOP here :: loop through all families here, plot dendrogram and PCA
# ===========
library(ggplot2)
library(qgraph)
ps <- ps.tweens.log   # use log normalized abandancies
#family.number <- twin.families[9]


for (family.number in twin.families[1:10]){
  print(family.number)
  twin.family.samples <- df.metadata.4timepoints[df.metadata.4timepoints$family_id==family.number, ]$file
  ps.onefamily <- phyloseq::subset_samples(ps, (sample_names(ps) %in% twin.family.samples))
  
  print(sample_data(ps.onefamily))
  
  # distance matrix
  ps.dist.w.unifrac <- phyloseq::distance(ps.onefamily, method="wunifrac", type="samples", fast=TRUE, parallel=TRUE)
  mat <- as(ps.dist.w.unifrac, "matrix")
  
  # create plot dendrogram
  plot(hclust(ps.dist.w.unifrac, method='ward.D2'))
  
  
  # plot PCA two components
  twin.ord <- phyloseq::ordinate(ps.onefamily, "NMDS", "unifrac")
  p2 <- phyloseq::plot_ordination(ps.onefamily, twin.ord, type="samples", color='twin_id', shape="human")
  print(p2)
  
  # plot grapg representation
  qgraph(ps.dist.w.unifrac, layout='spring', vsize=15)
  
}

## Add ggplot and combine all information on one page 
## check ordinate(ps.onefamily, "NMDS", "unifrac")



