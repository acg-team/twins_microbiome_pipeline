DATA=~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime/twin_demux-paired-end.qza
OUTPUT=~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime/demux.qzv
echo $DATA
echo $OUTPUT

qiime demux summarize --i-data $DATA --o-visualization $OUTPUT
qiime tools view $OUTPUT