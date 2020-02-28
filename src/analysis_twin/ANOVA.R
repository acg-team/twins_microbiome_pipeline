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

get_year <- function(date_factor){
  d <- as.Date(date_factor, format = "%Y-%m-%d")
  year(d)
}

# initialize matrix for intra / extra  distances
df.intr.extr<-data.frame("same_twin_diff_time","dif_twins_same_time","diff_twins_diff_time")


library(lubridate)

# CHANGE THE SOURCE PS FILI HERE 
ps <- ps.tweens.norm   # ps.tweens.norm  / ps.tweens.log
same.twin.different.years <- c()
same.year.different.twins <- c()
diff.year.different.twins <- c()

for (family.number in twin.families){
  print(family.number)
  twin.family.samples <- df.metadata.4timepoints[df.metadata.4timepoints$family_id==family.number, ]$file
  ps.onefamily <- phyloseq::subset_samples(ps, (sample_names(ps) %in% twin.family.samples))
  
  print(sample_data(ps.onefamily))
  
  # distance matrix
  ps.dist.w.unifrac <- phyloseq::distance(ps.onefamily, method="wunifrac", type="samples", parallel=TRUE)
  mat <- as(ps.dist.w.unifrac, "matrix")
  
  # create plot dendrogram
  #plot(hclust(ps.dist.w.unifrac, method='ward.D2'))
  
  # plot PCA two components
  twin.ord <- phyloseq::ordinate(ps.onefamily, method ="NMDS", distance ="unifrac")
  p2 <- phyloseq::plot_ordination(ps.onefamily, twin.ord, 
                                  type="samples", color='twin_id', label="collection_date",
                                  shape="human")
  
  # plot grapg representation
  p3 <- phyloseq::plot_net(ps.onefamily, ps.dist.w.unifrac, type = "samples",
                           point_label = "twin_id", maxdist = 1, 
                           color = "twin_id", point_size = 6)
  
  # add a table of samples
  tt <- ttheme_default(base_size = 6)
  
  tbl.meta <- tableGrob(df.metadata.4timepoints[df.metadata.4timepoints$family_id==family.number,], rows=NULL, theme=tt)
  tbl.dist <- tableGrob(mat, rows=NULL, theme=tt)
  
  g<-ggarrange(p2, p3, tbl.meta,  tbl.dist, 
            labels = c("ORD", "G"),
            ncol=2, nrow = 2,
            heights=c(2,2)
            )
  #print(g)
  ggsave(file=file.path(result_path, paste0(family.number,".jpg")), width=9, height = 5) 
  
  # for
  m <- df.metadata.4timepoints[df.metadata.4timepoints$family_id==family.number,]
  m["year"] <- lapply(m["collection_date"], get_year)
  twins <- unique(m$twin_id)
  years <- unique(m$year)
  
  
  # same twin different years
  for(twin in twins){
      sample <- m[(m$twin_id==twin),]$file
      
      # TODO: if we have 3 samples - allow it to be
      if (length(sample) < 2){  # if twin is only one - next
        next
      }
      
      print(paste(sample))
      d1 <- mat[ as.character(sample[1]), as.character(sample[2])]
      print(d1)
      same.twin.different.years <- c(same.twin.different.years,d1)
  }  
  
  # same year different twins
  for(year in years){
    sample <- m[(m$year==year),]$file
    if (length(sample) < 2){  # if twin is only one - next
      next
    }
    print(paste(sample))
    d2 <- mat[ as.character(sample[1]), as.character(sample[2])]
    print(d2)
    same.year.different.twins <- c(same.year.different.twins,d2)
  } 
  
  # dirrefent twins different years
  for( twin in twins){
    for (year in years){
      current_sample <- m[(m$twin_id == twin) & (m$year == year),]$file
      sample <- m[(m$twin_id != twin) & (m$year != year),]$file
      d3 <- mat[ as.character(current_sample), as.character(sample)]
      print(paste(current_sample, sample, d3))
      diff.year.different.twins <-c(diff.year.different.twins,d3)
    }
  }
  
  
}# end main loop
same.twin.different.years <- na.omit(same.twin.different.years)
same.year.different.twins <- na.omit(same.year.different.twins)
diff.year.different.twins <- na.omit(diff.year.different.twins)

print(paste("Time variability - Same twin different years (mean/std)", mean(same.twin.different.years), sd(same.twin.different.years)))
print(paste("Twins variability - Same year different twins (mean/std)", mean(same.year.different.twins), sd(same.year.different.twins)))
print(paste("Cross - Different year different twins (mean/std)", mean(diff.year.different.twins), sd(diff.year.different.twins)))



## check ordinate(ps.onefamily, "NMDS", "unifrac")
## filter out if only one example

# http://joey711.github.io/phyloseq/plot_ordination-examples
# https://github.com/joey711/phyloseq/wiki/ordinate








