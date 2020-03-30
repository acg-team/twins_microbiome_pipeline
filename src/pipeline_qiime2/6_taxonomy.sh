## Taxonomic analysis using greengenes
P=~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime/
GG=~/Projects_R/twins_microbiome_pipeline/tax_green_genes/
OUT=${P}taxonomy/



### First download the greengenes reference database on QIIME2 ###(https://docs.qiime2.org/2019.10/data-resources/)
### and import the otus file as a qiime2 object
# https://docs.qiime2.org/2020.2/data-resources/
qiime tools import \
    --type "FeatureData[Sequence]" \
    --input-path ${GG}gg_13_8_otus/rep_set/99_otus.fasta \
    --output-path ${OUT}99_otus.qza

### Then import the taxonomy file as a qiime2 object
qiime tools import --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path ${GG}gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
    --output-path ${OUT}ref-99taxonomy.qza


# TODO Do we need it?    
### Extracted the V4-V5 region by specifying the primers
#qiime feature-classifier 
#    extract-reads --i-sequences 99_otus.qza 
#    --p-f-primer CCAGCAGCYGCGGTAAN 
#    --p-r-primer CCGTCAATTNNTTTNANT 
#    --o-reads ref-99otus.qza
    

### Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ${OUT}99_otus.qza \
  --i-reference-taxonomy ${OUT}ref-99taxonomy.qza \
  --o-classifier ${OUT}classifier.qza 


## Run the classifier 
qiime feature-classifier classify-sklearn \
  --i-classifier ${OUT}classifier.qza \
  --i-reads ${P}otu/rep-seqs-dn-99.qza \
  --o-classification ${OUT}taxonomy.qza 

#qiime metadata tabulate \
#  --m-input-file taxonomy.qza \
#  --o-visualization taxonomy.qzv

