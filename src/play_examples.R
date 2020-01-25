library("phyloseq")
library("ggplot2")
library("doParallel")
library("foreach")
set.seed(711)

data(esophagus)
UniFrac(esophagus, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)

taxa_are_rows(esophagus)
taxa_are_rows(otu_table(esophagus))
