#################################
# TODO: add a nextflow pipeline?
#################################

source("src/load.R")
source("src/configure.R")
library(tidyverse)
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

load(file=file.path(metadata_path, metadata.file)) 

## Importing data
#This section describe how to import the data when they are not formatted with Casanova and chapcha...

#qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest --input-format PairedEndFastqManifestPhred33V2 --output-path demux-paired-end.qza

# we import data from filtered folder to reuse this dada2 functionality
# you need to import your data into QIIME 2 manually by first creating a “manifest file” 
#  and then using the qiime tools import command
# https://docs.qiime2.org/2020.2/tutorials/importing/

# 1 - create CSV manifest file
# sample-id     forward-absolute-filepath       reverse-absolute-filepath
qiime.manifest <- data.frame(sample.names, fnFs, fnRs)
colnames(qiime.manifest) <- c("sample-id", "forward-absolute-filepath", "reverse-absolute-filepath")
write_tsv(qiime.manifest, file.path(qiime_path, "manifest.tsv"))
print(file.path(qiime_path, "manifest.tsv"))
print(file.path(qiime_path, "twin_demux-paired-end.qza"))

# 2 - run import script 
# chmod +x 1_import.sh 2_vizualization.sh




