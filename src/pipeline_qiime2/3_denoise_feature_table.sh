P=~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime/

DATA=${P}twin_demux-paired-end.qza

OUTPUT_REP=${P}feature_table/rep-seqs-dada2.qza
OUTPUT_TABLE=${P}feature_table/table-dada2.qza
OUTPUT_STATS=${P}feature_table/stats-dada2.qza
VIZ_STATS=${P}feature_table/stats-dada2.qzv
EXPORT_DIR=${P}feature_table/rep-seqs-dada2


# repseqs.qza contains sequences of representative sequences identified by DADA2 and 
# q25table.qza is a table containing counts of each representative sequence in each sample.
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs $DATA \
    --p-trim-left-f 3 --p-trim-left-r 3 \
    --p-trunc-len-f 247 --p-trunc-len-r 235 \
    --o-representative-sequences $OUTPUT_REP \
    --o-table $OUTPUT_TABLE \
    --o-denoising-stats $OUTPUT_STATS
    
qiime metadata tabulate \
    --m-input-file $OUTPUT_STATS \
    --o-visualization $VIZ_STATS
    
qiime tools view $VIZ_STATS

qiime tools export --input-path $OUTPUT_REP --output-path $EXPORT_DIR
