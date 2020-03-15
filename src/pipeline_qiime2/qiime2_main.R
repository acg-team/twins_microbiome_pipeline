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

# 1 - create and save CSV manifest file
# sample-id     forward-absolute-filepath       reverse-absolute-filepath
#########################################################################
qiime.manifest <- data.frame(sample.names, fnFs, fnRs)
colnames(qiime.manifest) <- c("sample-id", "forward-absolute-filepath", "reverse-absolute-filepath")
write_tsv(qiime.manifest, file.path(qiime_path, "manifest.tsv"))

print(file.path(qiime_path, "manifest.tsv"))
print(file.path(qiime_path, "twin_demux-paired-end.qza"))

# 2 - run import script 
# chmod +x 1_import.sh 2_vizualization.sh 3_denoise_feature_table.sh


# 3 - Denoising and feture table construction


# 4 - OTU clustering

### get data into R
otus.99 <- qiime2R::read_qza("~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime/table-dn-99.qza")
otus.99$data
# add sequence names?

# 5 - tree

# 6 - Taxonomy with greengenes (in dada2 is with SILVA) - is that a difference?


### Build a phyloseq object
library(phyloseq)
tree<-read_qza("~/QIIME2/mvpics/rooted-tree.qza")

taxonomy<-read_qza("~/QIIME2/mvpics/taxonomy.qza")
tax_table<-do.call(rbind, strsplit(as.character(taxonomy$data$Taxon), "; "))
colnames(tax_table)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(tax_table)<-taxonomy$data$Feature.ID

metadata<-read.table("~/QIIME2/mvpics/sample-metadata.tsv", sep='\t', header=T, row.names=1, comment="")
metadata<-metadata[-1,]#remove the second line that specifies the data type

physeq<-phyloseq(otu_table(otus$data, taxa_are_rows = T), phy_tree(tree$data), tax_table(tax_table), sample_data(metadata))





