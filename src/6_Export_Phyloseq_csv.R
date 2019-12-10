
#### init: load packages and set path
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)
source("src/load_initialize.R")

### LOAD PREVIOUS DATA
load(file=file.path(files_intermediate, phyloseq.file)) 




##################### Export to Famuly_abundance.csv ############## 
# melt all sequences of one taxa to one abundance
# https://github.com/joey711/phyloseq/issues/418
glom <- phyloseq::tax_glom(ps, taxrank = 'Family')   # aglomerate taxa
dat <- psmelt(glom)   # convert to df
dat$Family <- as.character(dat$Family)  # convert Phylum to a character vector from a factor

Phylum_abundance <- aggregate(Abundance~Sample+Phylum, dat, FUN=sum)  # aggregate

# reorganize the table so that each phylum is a column
Phylum_abundance <- cast(Phylum_abundance, Sample ~ Phylum)

save(Phylum_abundance, file=file.path(result_path, "Phylum_abundance.RData")) 

write.table(Phylum_abundance, file.path(result_path,"Famuly_abundance.csv"), sep=",")