### use Phylogeq oblest (ps.tweens) 

#### init: load packages and set path
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)
source("src/load_initialize.R")

### LOAD PREVIOUS DATA
load(file=file.path(files_intermediate, phyloseq.file)) 



##############  EXPLORATORY ANALYSIS of Phyloseq object  #######
# TODO: show distribution of abandancies, max, min

# plot abundance of each taxa in all samples
plot_bar(ps.tweens, fill="Class")

# richness (number of taxa in each sample)
plot_richness(ps.tweens, x="BODY_SITE", color="Description")

# plot phylo tree with abandance
plot_tree(ps.tweens, color="Genus", size="abundance")

#####  Prevalence filtering
# Define prevalence of each taxa (in how many samples did each taxa appear at least once)
prev0 = apply(X = otu_table(ps.tweens),
              MARGIN = ifelse(taxa_are_rows(ps.tweens), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)}
              )
prevdf = data.frame(Prevalence = prev0, TotalAbundance = taxa_sums(ps.tweens), tax_table(ps.tweens))
View(head(prevdf))

# we can filter out Phyla with abundance less then 50
keepPhyla.gt50 = table(prevdf$Phylum)[(table(prevdf$Phylum) > 50)]
prevdf.Phyla.gt50 = subset(prevdf, Phylum %in% names(keepPhyla.gt50))
nrow(prevdf.Phyla.gt50)
View(head(prevdf.Phyla.gt50))

# we can define other rules for filering too

# we can filter out 5% of taxa with lowest abandances  
prevalenceThreshold = 0.05 * nsamples(ps.tweens) 
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
ps.tweens.1 = prune_taxa((prev0 > prevalenceThreshold), ps.tweens) 
ps.tweens.1

# Filter entries with unidentified Phylum.
ps.tweens.2 = subset_taxa(ps.tweens.1, Phylum %in% names(keepPhyla.gt50))
ps.tweens.2

# plot abandance of each Phyla
ggplot(prevdf.Phyla.gt50, aes(TotalAbundance, Prevalence, color = Phylum)) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 1, alpha = 0.7) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total~Abundance") +
  facet_wrap(~Phylum)





############ UNIFRAC Analysis   ##########
# UNIFRAC distance is calculated between pairs of samples (each sample represents an organismal community)
# - weighted (quantitative, abundance) and unweighted (qualitative, presence or absence)  


#########  Start exploration and analysis
# https://joey711.github.io/phyloseq-demo/unifrac.html
# do we need a rooted tree?

# Unifrac assess a distance btw two sets of microbial community based on their tree position and abandunce 
##### calculate UNIFRAC distance matric
# Takes a phyloseq object and method option,  and returns a distance object suitable for certain 
# ordination methods and other distance-based analyses
# methods: "unifrac": unweighted / UniFrac1, Unifrac2
# type: pairwise comparisons by samples

# 7 hours
tic()
ps.dist.unifrac <- phyloseq::distance(ps.tweens, method="unifrac", type="samples", fast=TRUE, parallel=TRUE)
toc()
save(ps.dist.unifrac, file=file.path(files_intermediate, phyloseq_analysis.file)) 


#### MDS analysis: heatMap
plot_heatmap(ps.tweens, distance = "unifrac", method="NMDS", sample.label="SampleType", taxa.label="Family")


plot(hclust(ps.dist.unifrac, method='ward.D'))


# do ordination (NMDS)
ord  <- ordinate(ps.tweens, "MDS", distance=ps.dist.unifrac)
plot(ord)









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