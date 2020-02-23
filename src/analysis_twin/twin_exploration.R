########### use combined Phyloseq object (ps.tweens) to analise microbial community

#### init: load packages and set path
source("src/configure.R")
source("src/load.R")
setwd(project_path)
print(conf)  # check current dataset


### LOAD PREVIOUS DATA
load(file=file.path(metadata_path, metadata.file))
load(file=file.path(files_intermediate, phyloseq.file))
# ps.tweens must be in the workspace now

# abanduncy normalization and log transform
ps.tweens.norm <- transform_sample_counts(ps.tweens, function(x) x / sum(x) )
ps.tweens.log <- transform_sample_counts(ps.tweens.norm, function(x) log(1 + x))

twin.families <- unique(df.metadata.4timepoints$family_id)

# FUNCTION: get a ps object for only one family
subset_by_family <- function(family_id){
  twin.family.samples <- df.metadata.4timepoints[df.metadata.4timepoints$family_id==twin.families[family.number], ]$file
  ps.onefamily <- subset_samples(ps.tweens.log, sample_names(ps.tweens) %in% twin.family.samples)
  ps.onefamily
}

print_family <- function(family_id){
  print(df.metadata.4timepoints[df.metadata.4timepoints$family_id==twin.families[family.number], ])
}

plot_bar_rich_tree <- function(ps.onefamily){
  # https://joey711.github.io/phyloseq/plot_bar-examples.html
  # 1 - bar plot -  ABUNDANCE of each taxa in all samples
  ggp2bar.family <- plot_bar(ps.onefamily, fill="Family")  
  ggsave(file=file.path(result_path, "bar_taxa_in_samples_Genus.png"), plot = ggp2bar.family, dpi = 300, width =30, height = 20)
  # todo - add ID to name or legend
  
  # 2 - richness - number of taxa in each sample (alpha-diversity)
  # TODO: study measures: observed, shennon etc
  ggp2rich <- plot_richness(ps.onefamily)   #x='family_id', measures="Observed"
  ggsave(file=file.path(result_path, "richness.png"), plot = ggp2rich, dpi = 300, width = 10, height = 8)
  
  
  # 3 - plot phylo tree with abandance
  # Incorrect number of arguments (7), expecting 5 for 'node_depth_edgelength'
  ggp2tree <- plot_tree(ps.onefamily, method = "treeonly", ladderize = "left", title = "Tree of twins taxa")
  ggsave(file=file.path(result_path, "tree.png"), plot = ggp2tree, dpi = 300, width = 20, height = 20)
  
  ggp2tree.aband <- plot_tree(ps.onefamily, color="Genus", size="abundance") 
  ggsave(file=file.path(result_path, "tree_abund.png"), plot = ggp2tree.aband, dpi = 300, width = 30, height = 30)
  
}


##############  EXPLORATORY ANALYSIS of Phyloseq object  #######
# ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
################################################################
# OTU = species
# ASV - an OTU with 100% similariry inside
# beta-diversity - type and quantity of OTU observed btw samples (sites) - what we are looking for
# alpha-diversity - diversity inside each sample (site)

# we are interested in beta-diversity, i.e. comparison of microbial composition of two samples, it is achived with
# - calculation of pairwise distances (distance matrix)
# - dimensional reduction (ordination) methods

# Unifrac assess a distance btw two sets of microbial community based on their tree position and abandunce 
# - weighted (quantitative, abundance) and unweighted (qualitative, presence or absence)  
# https://joey711.github.io/phyloseq-demo/unifrac.html


# LOOP::
family.number = 4
print_family(family.number)

ps.onefamily <- subset_by_family(family.number)

### calculate Unifrac DISTANCE matrix
# distances: "unifrac": unweighted / wunifrac-weighted
# type: pairwise comparisons by samples
# NOTE: ps.dist.unifrac - is a "dist" class (which package?) sutable for standard clustering analysis in R (hclust)
ps.family.dist.unifrac <- phyloseq::distance(ps.onefamily, method="wunifrac", type="samples", fast=TRUE, parallel=TRUE)

# HClustering linkage criteria or aglomerative methoe: "average" (= UPGMA), "centroid" (= UPGMC), "ward.D", "ward.D2"
plot(hclust(ps.family.dist.unifrac, method='ward.D'))

plot_heatmap(ps.onefamily, distance = "wunifrac", method="NMDS", sample.label="SampleType", taxa.label="Family")





# do ordination/dimensionality reduction (NMDS)
#ord  <- ordinate(ps.tweens, "MDS", distance=ps.dist.unifrac)
twin.ord <- ordinate(ps.onefamily, "NMDS", "unifrac")
p1 = plot_ordination(ps.onefamily, twin.ord, type="samples", color='twin_id')
print(p1)



twin.ord = ordinate(ps.onefamily, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps.onefamily, twin.ord, color="twin_id", shape="collection_date")

# TODO
# show distances 
# loop
# weited unifrac
# get a distance , calculated everage distance average 
# distance between same twin samples versus different twin samples
# run another fasta from quimm2








#########  1 - Preprocessing: filtering samples/taxa
# to distinguish twins we have to filter out 
# - taxa with small variance (which are stable) - common
# - shall it be Phyla? Class? which one?
# - filter out stable taxa?


# Prevalence filtering: Define prevalence of each taxa (in how many samples did each taxa appear at least once)
# Fast: all in 5 min
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
#note: we can define other rules for filering too

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
ggp2.prevalence <- ggplot(prevdf.Phyla.gt50, aes(TotalAbundance, Prevalence, color = Phylum)) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 1, alpha = 0.7) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total~Abundance") +
  facet_wrap(~Phylum)
ggsave(file=file.path(result_path, "prevalence.pdf"), plot = ggp2.prevalence, dpi = 300, width = 30, height = 20)



