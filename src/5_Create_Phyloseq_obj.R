## Combining a S4 Phyloseq object for further manipulation
## phyloseq object is an experiment level data structure
## http://rpubs.com/lgschaerer/515637
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


############    BUILD a Phyloseq object
# Phyloseq obgect is a  typical amplicon sequencing experiment in one single data object 

# we need all sampleID to be the same in rownames(seqtab) and samples metadata (df.metadata$SampleID)
# we need a SampleID in order to Phyloseq object be valid (change the name of /file/ column to SampleID)
names(df.metadata)[names(df.metadata)=="file"] <- "SampleID"
# assign the names of samples (ERR138...) to metadata rows instead of 1,2,3...
rownames(df.metadata) <- df.metadata$SampleID

# should be all TRUE (sanity check)
all(rownames(seqtab) %in% df.metadata$SampleID) 

# build phyloseq object
feature.table <- otu_table(seqtab, taxa_are_rows = FALSE)
metadata.table <- sample_data(df.metadata)
tree.final <- phy_tree(treeNJ)   # fitGTR$tree # phylo object (fitGTR$tree) we temporarily use NJ tree here instead of RAXML: 

# assign long sequence names back (we made it short after RAxML)
taxa_names(tree.final) <- colnames(feature.table)

# control
rownames(metadata.table)
colnames(feature.table)
taxa_names(tree.final)

# Create an object 
ps.tweens <- phyloseq::phyloseq(
              tax_table(taxtab), 
              metadata.table,
              feature.table,
              tree.final
              )

save(ps.tweens, file=file.path(files_intermediate, phyloseq.file)) 

## CHECK UP: phyloseq-class experiment-level object
otu_table(ps)     # OTU Table:         [ 8299 taxa and 3288 samples ]
sample_data(ps)   # Sample Data:       [ 3288 samples by 8 features ]
tax_table(ps)     # Taxonomy Table:    [ 8299 taxa by 6 taxonomic ranks ]
phy_tree(ps)      # Phylogenetic Tree: [ 8299 tips, ??? internal nodes ] - need to re run the whole workflow








