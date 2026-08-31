#!/usr/bin/env bash
# =============================================================================
# Core biom
# =============================================================================

set -Eeuo pipefail
shopt -s nullglob

# ==================== CONFIGURATION ====================

ROOTDIR="/nvme/bio/data_fungi/valormicro_V2"
PROJECT_DIR="${ROOTDIR}"
RESULTS_DIR="${PROJECT_DIR}/02_amplicon_pipeline"
DATABASE=/nvme/bio/data_fungi/valormicro_V2/02_amplicon_pipeline/04_database_files

NTHREADS=16
TMPDIR="${ROOTDIR}/tmp"
QIIME_ENV="qiime2-amplicon-2025.7"

DATABASE_DIR="${ROOTDIR}/98_databasefiles"
QIIME_CORE_DIR="${RESULTS_DIR}/05_qiime2/core"
LOG_DIR="${RESULTS_DIR}/logs"

REP_SEQS="${QIIME_CORE_DIR}/rep-seqs.qza"
TAXONOMY="${QIIME_CORE_DIR}/taxonomy.qza"

mkdir -p "${TMPDIR}" "${DATABASE_DIR}" "${LOG_DIR}"
export TMPDIR

# ==================== FONCTIONS ====================

log() {
    printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOG_DIR}/03_taxonomy.log"
}

die() {
    log "ERREUR : $*"
    exit 1
}

trap 'rc=$?; log "ERREUR : code ${rc}, ligne ${LINENO} : ${BASH_COMMAND}"; exit "${rc}"' ERR

# ==================== CONTROLES ====================

log "Debut du core biom"
log "Environnement QIIME2 : ${QIIME_ENV}"
log "Nombre de jobs pour classify-sklearn : ${NTHREADS}"

command -v conda >/dev/null 2>&1 || die "Conda est introuvable dans le PATH."

conda run -n "${QIIME_ENV}" qiime --version \
    || die "QIIME2 ne demarre pas dans l'environnement ${QIIME_ENV}."


cd "${RESULTS_DIR}/05_qiime2"

conda run -n "${QIIME_ENV}" qiime diversity core-metrics-phylogenetic \
       --i-phylogeny tree/rooted-tree.qza \
       --i-table core/table.qza \
       --p-sampling-depth 3547 \
       --m-metadata-file $DATABASE/sample-metadata.tsv \
       --o-rarefied-table core/RarTable.qza \
       --o-observed-features-vector core/Vector-observed_asv.qza \
       --o-shannon-vector core/Vector-shannon.qza \
       --o-evenness-vector core/Vector-evenness.qza \
       --o-faith-pd-vector core/Vector-faith_pd.qza
     
conda run -n "${QIIME_ENV}" qiime feature-table core-features \
        --i-table core/RarTable.qza \
        --p-min-fraction 0.1 \
        --p-max-fraction 1.0 \
        --p-steps 10 \
        --o-visualization visual/CoreBiom-all.qzv  
        
conda run -n "${QIIME_ENV}" qiime tools export --input-path core/RarTable.qza --output-path export/core/RarTable   
conda run -n "${QIIME_ENV}" qiime tools export --input-path visual/CoreBiom-all.qzv --output-path export/visual/CoreBiom-all
biom convert -i export/core/RarTable/feature-table.biom -o export/core/RarTable/table-from-biom.tsv --to-tsv
sed '1d ; s/\#OTU ID/ASV_ID/' export/core/RarTable/table-from-biom.tsv > export/core/RarTable/ASV.tsv

conda run -n "${QIIME_ENV}" qiime tools export --input-path core/Vector-faith_pd.qza --output-path export/core/Vector-faith_pd
conda run -n "${QIIME_ENV}" qiime tools export --input-path core/Vector-evenness.qza --output-path export/core/Vector-evenness
conda run -n "${QIIME_ENV}" qiime tools export --input-path core/Vector-shannon.qza --output-path export/core/Vector-shannon
conda run -n "${QIIME_ENV}" qiime tools export --input-path core/Vector-observed_asv.qza --output-path export/core/VVector-observed_asv

conda run -n "${QIIME_ENV}" qiime tools export --input-path core/taxonomy.qza --output-path export/taxonomy
