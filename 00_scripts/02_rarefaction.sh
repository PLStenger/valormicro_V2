#!/usr/bin/env bash

WORKING_DIRECTORY=/nvme/bio/data_fungi/valormicro_V2/02_amplicon_pipeline/05_qiime2
DATABASE=/nvme/bio/data_fungi/valormicro_V2/02_amplicon_pipeline/04_database_files

# Aim: rarefy a feature table to compare alpha/beta diversity results
# A good forum to understand what it does :
# https://forum.qiime2.org/t/can-someone-help-in-alpha-rarefaction-plotting-depths/4580/16

cd $WORKING_DIRECTORY

eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2025.7

qiime diversity alpha-rarefaction \
--i-table core/table_min2samples.qza \
--i-phylogeny tree/rooted-tree.qza \
  --p-max-depth 12618 \
  --p-min-depth 1 \
  --m-metadata-file $DATABASE/sample-metadata.tsv \
  --o-visualization visual/alpha-rarefaction.qzv
