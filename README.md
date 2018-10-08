# twins_microbiome_pipeline
Bioinformatics procession of raw Illumina reads for microbiome analysis of UK twins study with dada2 pipeline

Contains a complete workflow from extracting metainformation to assigning taxa information and phylogenetic tree of beta diversion.
Pipeline as from https://www.bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html


## raw Data and third-party data
- download all raw illumina reads from http://www.ebi.ac.uk/ena/data/view/PRJEB13747
- place them into /data/raw

## run subsequently
- 1_metadata.R - that script automatically downloads all metainformation atatched to every sample, age, family_id etc. If nessesary, modify the script to get more features then currently present

- 2_BIGD_dada_SV_table.R - the script implementing dada2 pipeline to process raw reads into SV table (analogius to OTU, but on variant level, not species). That is very long running script (~ 40 hours on a powerful laptop)