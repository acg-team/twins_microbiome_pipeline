#################################
# TODO: add a nextflow pipeline?
#################################

source("src/load.R")
source("src/configure.R")
library(tidyverse)
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

load(file=file.path(metadata_path, metadata.file)) 


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


##### Here we have to run qiime worklow by subsequently running 6 scripts
# chmod +x 1_import.sh 2_vizualization.sh 3_denoise_feature_table.sh
# 2 - run import script 
# 3 - Denoising and feture table construction
# 4 - OTU clustering
# 5 - tree
# 6 - Taxonomy with greengenes (in dada2 is with SILVA) - is that a difference?



############ Build a phyloseq object
library(phyloseq)
qiime.data.path <- "~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime"

# form meatadata
qiime.metadata <- df.metadata
names(qiime.metadata)[names(qiime.metadata)=="file"] <- "SampleID"

# form feature table
otus.99 <- qiime2R::read_qza(file.path(qiime.data.path, 'otu/table-dn-99.qza'))
feature.table <- otus.99$data

# tree
tree <- qiime2R::read_qza(file.path(qiime.data.path, 'tree/rooted-tree.qza'))
tree$data

# taxonomy --- TODO - something wrong here. we have to run slassifier on data?
tax.table <- qiime2R::read_qza(file.path(qiime.data.path, 'taxonomy/classifier.qza'))

# TODO here - !

taxonomy<-read_qza("~/QIIME2/mvpics/taxonomy.qza")
tax_table<-do.call(rbind, strsplit(as.character(taxonomy$data$Taxon), "; "))
colnames(tax_table)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(tax_table)<-taxonomy$data$Feature.ID




physeq<-phyloseq(
  otu_table(feature.table, taxa_are_rows = T), ### is that correct?
  phy_tree(tree$data), 
  tax_table(tax_table), 
  sample_data(metadata)
  )

# save it




