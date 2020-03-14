P=~/Projects_R/twins_microbiome_pipeline/data_set_twin/raw/qiime/
REP=${P}feature_table/rep-seqs-dada2.qza

OUT_ALIGNMENT=${P}tree/aligned-rep-seqs.qza
OUT_MASKED_ALIGN=${P}tree/masked-aligned-rep-seqs.qza
OUT_UNROOTED_TREE=${P}tree/unrooted-tree.qza
OUT_ROOTED_TREE=${P}tree/rooted-tree.qza


## Generate a tree for phylogenetic diversity analyses

qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences $REP \
    --o-alignment $OUT_ALIGNMENT \
    --o-masked-alignment $OUT_MASKED_ALIGN \
    --o-tree $OUT_UNROOTED_TREE \
    --o-rooted-tree $OUT_ROOTED_TREE