P=~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime/

TABLE=${P}feature_table/table-dada2.qza
REP=${P}feature_table/rep-seqs-dada2.qza

OUTPUT_CLUST_TABLE99=${P}otu/table-dn-99.qza
OUTPUT_CLUST_SEQ99=${P}otu/rep-seqs-dn-99.qza
EXPORT_DIR99=${P}otu/99otu_table

OUTPUT_CLUST_TABLE97=${P}otu/table-dn-97.qza
OUTPUT_CLUST_SEQ97=${P}otu/rep-seqs-dn-97.qza
EXPORT_DIR97=${P}otu/97otu_table

# De novo clustering using 99% threshold  
qiime vsearch cluster-features-de-novo \
    --i-table $TABLE \
    --i-sequences $REP \
    --p-perc-identity 0.99 \
    --o-clustered-table $OUTPUT_CLUST_TABLE99 \
    --o-clustered-sequences $OUTPUT_CLUST_SEQ99
    
# Export the feature table
qiime tools export --input-path $OUTPUT_CLUST_TABLE99 --output-path $EXPORT_DIR99


## (Optional) de novo OTU clustering - in the example below, using 97% threshold  
qiime vsearch cluster-features-de-novo \
    --i-table $TABLE \
    --i-sequences $REP \
    --p-perc-identity 0.97 \
    --o-clustered-table $OUTPUT_CLUST_TABLE97 \
    --o-clustered-sequences $OUTPUT_CLUST_SEQ97
    
## Export the feature table
qiime tools export --input-path $OUTPUT_CLUST_TABLE97 --output-path $EXPORT_DIR97
    
    