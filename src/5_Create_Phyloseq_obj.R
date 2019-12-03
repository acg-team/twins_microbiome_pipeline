## Combining a S4 Phyloseq object for further manipulation
## phyloseq object is an experiment level data structure
##########################################################################

#### init: load packages and set path
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)
source("src/load_initialize.R")

### LOAD PREVIOUS DATA
load(file=file.path(metadata_path, metadata.file))
load(file=file.path(files_intermediate, seqtab.file)) 
load(file=file.path(files_intermediate, seqtab.snames.file)) 
load(file=file.path(files_intermediate, taxtab.file))
load(file=file.path(files_intermediate, phylo.file)) 


#### BUILD a Phyloseq obgect which carries all information about tree in one file

# change the name of /file/ column to SampleID
names(df.metadata)[names(df.metadata)=="file"] <- "SampleID"
# assign the names of samples (ERR138...) to rows instead of 1,2,3...
rownames(df.metadata) <- df.metadata$SampleID
# should be all TRUE (sanity check)
all(rownames(seqtab) %in% df.metadata$SampleID) 

# build phyloseq object
ps <- phyloseq::phyloseq(
              tax_table(taxtab), 
              sample_data(df.metadata),
              otu_table(seqtab, taxa_are_rows = FALSE), 
              phy_tree(treeNJ)   # phylo object (fitGTR$tree) we temporarily use NJ tree here instead of RAXML: 
              )

save(ps, file=file.path(files_intermediate, phyloseq.file)) 


##### Phyloseq TODO here




#####################
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

# TODO - put twins together


