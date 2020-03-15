## Taxonomic analysis using greengenes
P=~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime/
GG=~/Projects_R/twins_microbiome_pipeline/green_genes/
OUT=${P}taxonomy/



### First download the greengenes reference database on QIIME2 ###(https://docs.qiime2.org/2019.10/data-resources/)
### and import the otus file as a qiime2 object
# https://docs.qiime2.org/2020.2/data-resources/
qiime tools import 
    --type "FeatureData[Sequence]" 
    --input-path ${GG}gg_13_8_otus/rep_set/99_otus.fasta 
    --output-path ${OUT}99_otus.qza

### Then import the taxonomy file as a qiime2 object
qiime tools import --type 'FeatureData[Taxonomy]' 
    --input-format HeaderlessTSVTaxonomyFormat 
    --input-path ${GG}gg_13_8_otus/taxonomy/99_otu_taxonomy.txt 
    --output-path ${OUT}ref-99taxonomy.qza


# TODO from here?    
### Extracted the V4-V5 region by specifying the primers
qiime feature-classifier 
    extract-reads --i-sequences 99_otus.qza 
    --p-f-primer CCAGCAGCYGCGGTAAN 
    --p-r-primer CCGTCAATTNNTTTNANT 
    --o-reads ref-99otus.qza
    
### Train the classifier

    qiime feature-classifier fit-classifier-naive-bayes
  --i-reference-reads ref-99otus.qza 
  --i-reference-taxonomy ref-99taxonomy.qza 
  --o-classifier classifier.qza

## Test the classifier
We verify that the classifier works by classifying the representative sequences from Dobay et al.

qiime feature-classifier classify-sklearn 
  --i-classifier classifier.qza 
  --i-reads rep-seqs-dn-99.qza
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv