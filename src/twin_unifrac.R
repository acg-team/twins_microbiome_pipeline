source("src/load.R")
source("src/configure.R")

load(file=file.path(metadata_path, metadata.file))
load(file=file.path(files_intermediate, phyloseq.file))
# ps.tweens must be in the workspace now

head(df.metadata.4timepoints)
ps.tweens

### calculate Unifrac DISTANCE matrix
# distances: "unifrac": unweighted / wunifrac-weighted
# type: pairwise comparisons by samples
# NOTE: ps.dist.unifrac - is a "dist" class (which package?) sutable for standard clustering analysis in R (hclust)

tic()
ps.dist.w.unifrac <- phyloseq::distance(ps.tweens, method="wunifrac", type="samples", fast=TRUE, parallel=TRUE)
toc()
save(ps.dist.w.unifrac, file=file.path(alalysis_path, "unifrac.RData")) 

tic()
ps.dist.unifrac <- phyloseq::distance(ps.tweens, method="unifrac", type="samples", fast=TRUE, parallel=TRUE)
toc()
save(ps.dist.w.unifrac,ps.dist.unifrac, file=file.path(alalysis_path, "unifrac.RData")) 

attributes(ps.dist.unifrac)
mat <- as(ps.dist.unifrac, "matrix")

### EXlopration
plot(hclust(ps.dist.w.unifrac, method='ward.D2'))
# You probably forgot to remove the zero variance columns in your matrix
# https://support.bioconductor.org/p/45944/



