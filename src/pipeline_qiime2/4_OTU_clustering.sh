DATA=~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime/twin_demux-paired-end.qza
OUTPUT_REP=~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime/rep-seqs-dada2.qza


# De novo clustering using 99% threshold  

qiime vsearch cluster-features-de-novo \
    --i-table table.qza \
    --i-sequences rep-seqs.qza \
    --p-perc-identity 0.99 \
    --o-clustered-table table-dn-99.qza \
    --o-clustered-sequences rep-seqs-dn-99.qza
    
# Export the feature table
qiime tools export   
    --input-path table-dn-99.qza 
    --output-path 99otu_table