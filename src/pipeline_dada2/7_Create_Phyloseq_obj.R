## Combining a S4 Phyloseq object for further manipulation
## phyloseq object is an experiment level data structure
## http://rpubs.com/lgschaerer/515637
## https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
##########################################################################


### LOAD PREVIOUS DATA
load(file=file.path(metadata_path, metadata.file))
load(file=file.path(files_intermediate_dada, seqtab.file)) 
load(file=file.path(files_intermediate_dada, seqtab.snames.file)) 
load(file=file.path(files_intermediate_dada, tax.fname))
load(file=file.path(files_intermediate_dada, phylo.file)) 
load(file=file.path(files_intermediate_dada, msa.file)) 



###### sanity check and format adjustment
  # we need all sampleID to be the same in rownames(seqtab) and samples metadata (df.metadata$SampleID)
  # we need a SampleID in order to Phyloseq object be valid (change the name of /file/ column to SampleID)

if(conf$dataset == "TWIN"){
  names(df.metadata)[names(df.metadata)=="file"] <- "SampleID"
} else if(conf$dataset == "BODYFL"){
  # already have SampleID as a column name
  head(df.metadata)
}

# assign the names of samples (ERR138...) to metadata rows instead of 1,2,3...
rownames(df.metadata) <- df.metadata$SampleID

# must be TRUE : sample ID shall be the same
if( !all(rownames(seqtab) %in% df.metadata$SampleID)){
  stop(" Sample names must coinside!")
}

# metadate shall be presented for every sample
for(seqname in rownames(seqtab)){
  if(seqname %in% df.metadata$SampleID){
  } else{
    stop(paste("missed metadata for :", seqname, " sample"))
  }
}

rowSums(seqtab) # check for non-zero




############    BUILD a Phyloseq object
# Phyloseq obgect is a  typical amplicon sequencing experiment in one single data object 
# canonical OTU mast be taxa x samples?  Here samples x taxa, so we set taxa_are_rows = FALSE
# https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

head(asv_sequences)  # here should be short_name -> asv sequences correcpondence
shortname_to_asv <- names(asv_sequences)
names(shortname_to_asv) <- asv_sequences
head(shortname_to_asv)

metadata.table <- sample_data(df.metadata)

# make sure the row.names are the same for asll tables (shall be otu/asv)
# we will use short name, like names(asv_sequences)
taxonomy.table <- phyloseq::tax_table(as.matrix(taxtab))
taxa_names(taxonomy.table)


feature.table <- otu_table(as.matrix(seqtab), taxa_are_rows = FALSE)  # seqtab = ERR128(row) x TCGA(cols, taxa)
# substitute the long names with the short names
namesasv <- taxa_names(feature.table)
new_names <- lapply(namesasv, function(x) {shortname_to_asv[x]})
taxa_names(feature.table) <- new_names
taxa_names(feature.table)


tree.final <-  my.tree$bestTree   # raxml is rooted / GTP is unrooted
taxa_names(tree.final)



# Create a ps object 
ps <- phyloseq::phyloseq(
              taxonomy.table, 
              metadata.table,
              feature.table,
              tree.final
              )


# Sanity: Check if there are ASVs with no counts
any(taxa_sums(ps) == 0)
is.rooted(phy_tree(ps))
any(is.na(ps@otu_table@.Data))


## CHECK UP: phyloseq-class experiment-level object
# https://joey711.github.io/phyloseq/import-data.html
dim( otu_table(ps) )   # OTU Table:      [ 8299 taxa and 3288 samples ] - You must also specify if the species are rows or columns
dim( sample_data(ps) ) # Sample Data:    [ 3288 samples by 8 features ] - rownames must match the sample names in the otu_table
dim( tax_table(ps) )   # Taxonomy Table: [ 8299 taxa by 6 taxonomic ranks ] - The rownames must match the OTU names (taxa_names) 
phy_tree(ps)     # Phylogenetic Tree: [ 8299 tips, ??? internal nodes ] - need to re run the whole workflow


# form a file name depending on tools and parameter used
folder.suffix <- paste0(
  conf$dataset, "_", conf$pipeline, 
  "_Q", dada_param$QUALITY_THRESHOLD, 
  "_mEE", dada_param$maxEE[1], dada_param$maxEE[2], 
  "_trL", dada_param$trimLeft[1], dada_param$trimLeft[2],
  "_trR", dada_param$trimRight[1], dada_param$trimRight[2],
  "_truncLn", dada_param$truncLen[1], "_", dada_param$truncLen[2],
  "_msa_", tools_param$MSA_aligner,
  "_tax_", tools_param$tax_method
)

file.suffix <- paste0(folder.suffix, ".RData")
phyloseq.file <- paste0("run_", file.suffix)



# filter.log shall come from 3 file
save(ps, filter.log, conf, dada_param, tools_param , my.msa, my.tree, file=file.path(files_intermediate_dada, phyloseq.file)) 





