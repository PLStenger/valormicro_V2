#!/usr/bin/env bash

#set -euo pipefail

# ==================== CONFIGURATION ====================
export ROOTDIR="/nvme/bio/data_fungi/valormicro_V2"
export NTHREADS=16
export TMPDIR="${ROOTDIR}/tmp"
export QIIME_ENV="qiime2-amplicon-2025.7"
export FASTQC_ENV="fastqc"
export TRIMMOMATIC_ENV="trimmomatic"

mkdir -p "$TMPDIR"


eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2025.7


PROJECT_DIR="/nvme/bio/data_fungi/valormicro_V2"
RESULTS_DIR="${PROJECT_DIR}/02_amplicon_pipeline"


# ==================== 07 CLASSIFICATION ====================
log "Classification taxonomique SILVA 138.2"

cd "${PROJECT_DIR}/98_databasefiles"

CLASSIFIER_URL="/nvme/bio/data_fungi/valormicro_nc/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
CLASSIFIER="${ROOTDIR}/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"

if [ ! -f "$CLASSIFIER" ]; then
    cp "$CLASSIFIER_URL" "$CLASSIFIER" || { log "ERREUR: Classifier manquant"; exit 1; }
fi

conda run -n "$QIIME_ENV" qiime tools validate "$CLASSIFIER" || { log "ERREUR: Classifier invalide"; exit 1; }


cd "${RESULTS_DIR}/05_qiime2/core"

conda run -n "$QIIME_ENV" qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza \
    --p-n-jobs 16

log "✅ Classification réussie"
