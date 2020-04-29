################################################
# QIIME analysys:: 
# TWIN dataset - ready
# BDF - TODO!
################################################

source("src/load.R")
source("src/configure.R")
library(tidyverse)
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)


print(metadata_path)
print(qiime_path)
load(file=file.path(metadata_path, metadata.file)) 
df.metadata


# you need to import your data into QIIME 2 manually by first creating a “manifest file” 
#  and then using the qiime tools import command
# https://docs.qiime2.org/2020.2/tutorials/importing/
# https://drive.google.com/file/d/1dtl1XpvvNO015C2rIMEsNWm0u6cjQO6-/view?ts=5e625bda

# 1 - create and save CSV manifest file
# sample-id     forward-absolute-filepath       reverse-absolute-filepath
#########################################################################
qiime.manifest <- data.frame(sample.names, fnFs, fnRs)
colnames(qiime.manifest) <- c("sample-id", "forward-absolute-filepath", "reverse-absolute-filepath")
write_tsv(qiime.manifest, file.path(qiime_path, "manifest.tsv"))

print(file.path(qiime_path, "manifest.tsv"))
print(file.path(qiime_path, "twin_demux-paired-end.qza"))


##### Here we have to run qiime worklow by subsequently running 6 bash scripts in shell
# chmod +x 1_import.sh 2_vizualization.sh 3_denoise_feature_table.sh
# 2 - run import script 
# 3 - Denoising and feture table construction
# 4 - OTU clustering
# 5 - tree
# 6 - Taxonomy with greengenes (in dada2 is with SILVA) - is that a difference?



############ Build a phyloseq object
library(phyloseq)
#qiime.data.path <- "~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime"

# form meatadata
qiime.metadata <- df.metadata
names(qiime.metadata)[names(qiime.metadata)=="file"] <- "SampleID"

# read a  feature table (OTU)
# for OTU in Phyloseq we have to check if taxa are rows and set flag taxa_are_rows = TRUE / FALSE
otus.99 <- qiime2R::read_qza(file.path(qiime_path, 'otu/table-dn-99.qza'))
feature.table.qiime.matr <- otus.99$data   # variants(TAXA) x samples(columns)  #here it is matrix
feature.table.qiime <- otu_table(feature.table.qiime.matr, taxa_are_rows = TRUE)

# tree
tree.qiime <- qiime2R::read_qza(file.path(qiime_path, 'tree/rooted-tree.qza'))
tree.qiime$data

# taxonomy 
# Taxonomy Table: [ 287 taxa by 7 taxonomic ranks ] ??
taxonomy.qiime <- qiime2R::read_qza(file.path(qiime_path, 'taxonomy/taxonomy.qza')) # TODO!

# form a taxonomy from qiime object
tax_table <- do.call(rbind, strsplit(as.character(taxonomy.qiime$data$Taxon), "; "))
colnames(tax_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(tax_table) <- taxonomy.qiime$data$Feature.ID


###### form a phyloseq object
ps.twins.qiime <- phyloseq(
  feature.table.qiime, # done
  phy_tree(tree.qiime$data), 
  tax_table(tax_table), # done
  sample_data(qiime.metadata) # done
  )

##################### save it here
save(ps.twins.qiime, file=file.path(files_intermediate, 'phyloseq_object_qiime.RData')) 

# SANITY CHECK
ps.twins.qiime
# OTU Table:         [ taxa by  samples ] - You must also specify if the species are rows or columns
# Sample Data:       [ samples by features ] - rownames must match the sample names in the otu_table
# Taxonomy Table:    [ taxa by taxonomic ranks ] - The rownames must match the OTU names (taxa_names) 
# Phylogenetic Tree: [ tips, ??? internal nodes ] - need to re run the whole workflow




