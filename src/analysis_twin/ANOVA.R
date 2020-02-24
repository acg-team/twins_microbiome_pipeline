# This is a fule for AKOS to modify
# https://github.com/joey711/phyloseq/wiki/ordinate

library(ggplot2)
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


ps <- ps.tweens.norm   # use log normalized abandancies
#family.number <- twin.families[1]


for (family.number in twin.families[1:50]){
  print(family.number)
  twin.family.samples <- df.metadata.4timepoints[df.metadata.4timepoints$family_id==family.number, ]$file
  ps.onefamily <- phyloseq::subset_samples(ps, (sample_names(ps) %in% twin.family.samples))
  
  print(sample_data(ps.onefamily))
  
  # distance matrix
  ps.dist.w.unifrac <- phyloseq::distance(ps.onefamily, method="wunifrac", type="samples", fast=TRUE, parallel=TRUE)
  mat <- as(ps.dist.w.unifrac, "matrix")
  
  # create plot dendrogram
  #plot(hclust(ps.dist.w.unifrac, method='ward.D2'))
  
  # plot PCA two components
  twin.ord <- phyloseq::ordinate(ps.onefamily, "NMDS", "unifrac")
  p2 <- phyloseq::plot_ordination(ps.onefamily, twin.ord, type="samples", color='twin_id', shape="human")
  
  # plot grapg representation
  p3 <- phyloseq::plot_net(ps.onefamily, point_label = "twin_id", maxdist = 1.5, color = "twin_id")
  
  # add a table of samples
  tt <- ttheme_default(base_size = 6)
  
  tbl <- tableGrob(df.metadata.4timepoints[df.metadata.4timepoints$family_id==family.number,], rows=NULL, theme=tt)
  
  g<-ggarrange(p2, p3, tbl,
            labels = c("ORD", "G"),
            ncol=2, nrow = 2,
            heights=c(2,1)
            )
  #print(g)
  ggsave(file=file.path(result_path, paste0(family.number,".jpg")), width=9, height = 4) 
}


## check ordinate(ps.onefamily, "NMDS", "unifrac")
## filter out if only one example








