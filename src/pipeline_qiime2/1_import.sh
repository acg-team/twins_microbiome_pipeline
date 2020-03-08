DATA=~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime/manifest.tsv
OUTPUT=~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime/twin_demux-paired-end.qza

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path $DATA --input-format PairedEndFastqManifestPhred33V2 --output-path $OUTPUT
