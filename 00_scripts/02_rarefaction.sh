#!/usr/bin/env bash

WORKING_DIRECTORY=/nvme/bio/data_fungi/valormicro_V2/02_amplicon_pipeline/05_qiime2
DATABASE=/nvme/bio/data_fungi/valormicro_V2/02_amplicon_pipeline/04_database_files

# Aim: rarefy a feature table to compare alpha/beta diversity results
# A good forum to understand what it does :
# https://forum.qiime2.org/t/can-someone-help-in-alpha-rarefaction-plotting-depths/4580/16

cd $WORKING_DIRECTORY

eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2025.7

##########################################################################################
# Arbre sans contingence

   qiime alignment mafft \
        --i-sequences core/rep-seqs.qza \
        --p-n-threads 8 \
        --o-alignment tree/aligned-rep-seqs.qza

    qiime alignment mask \
        --i-alignment tree/aligned-rep-seqs.qza \
        --o-masked-alignment tree/masked-aligned-rep-seqs.qza

    qiime phylogeny fasttree \
        --i-alignment tree/masked-aligned-rep-seqs.qza \
        --o-tree tree/unrooted-tree.qza

    qiime phylogeny midpoint-root \
        --i-tree tree/unrooted-tree.qza \
        --o-rooted-tree tree/rooted-tree.qza
##########################################################################################
        

qiime diversity alpha-rarefaction \
--i-table core/table.qza \
--i-phylogeny tree/rooted-tree.qza \
  --p-max-depth 15959 \
  --p-min-depth 1 \
  --m-metadata-file $DATABASE/sample-metadata.tsv \
  --o-visualization visual/alpha-rarefaction.qzv
