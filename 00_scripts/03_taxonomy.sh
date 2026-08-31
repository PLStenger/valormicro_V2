#!/usr/bin/env bash
# =============================================================================
# Classification taxonomique des ASV avec le classifieur SILVA 138.2
# =============================================================================

set -Eeuo pipefail
shopt -s nullglob

# ==================== CONFIGURATION ====================

ROOTDIR="/nvme/bio/data_fungi/valormicro_V2"
PROJECT_DIR="${ROOTDIR}"
RESULTS_DIR="${PROJECT_DIR}/02_amplicon_pipeline"

NTHREADS=16
TMPDIR="${ROOTDIR}/tmp"
QIIME_ENV="qiime2-amplicon-2025.7"

DATABASE_DIR="${ROOTDIR}/98_databasefiles"
QIIME_CORE_DIR="${RESULTS_DIR}/05_qiime2/core"
LOG_DIR="${RESULTS_DIR}/logs"

CLASSIFIER_SOURCE="/nvme/bio/data_fungi/valormicro_nc/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
CLASSIFIER="${DATABASE_DIR}/silva-138.2-ssu-nr99-515f-926r-classifier.qza"

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

log "Debut de la classification taxonomique SILVA 138.2"
log "Environnement QIIME2 : ${QIIME_ENV}"
log "Nombre de jobs pour classify-sklearn : ${NTHREADS}"

command -v conda >/dev/null 2>&1 || die "Conda est introuvable dans le PATH."

conda run -n "${QIIME_ENV}" qiime --version \
    || die "QIIME2 ne demarre pas dans l'environnement ${QIIME_ENV}."

[[ -f "${REP_SEQS}" ]] \
    || die "Sequences representatives absentes : ${REP_SEQS}"

# ==================== CLASSIFIEUR ====================

if [[ ! -f "${CLASSIFIER}" ]]; then
    log "Copie du classifieur SILVA vers ${DATABASE_DIR}"

    [[ -f "${CLASSIFIER_SOURCE}" ]] \
        || die "Classifieur source introuvable : ${CLASSIFIER_SOURCE}"

    cp -v "${CLASSIFIER_SOURCE}" "${CLASSIFIER}"
fi

[[ -f "${CLASSIFIER}" ]] \
    || die "Classifieur introuvable apres copie : ${CLASSIFIER}"

log "Validation de l'artefact classifieur"

conda run -n "${QIIME_ENV}" qiime tools validate "${CLASSIFIER}" \
    || die "Classifieur QIIME2 invalide : ${CLASSIFIER}"

# ==================== CLASSIFICATION ====================

if [[ -f "${TAXONOMY}" ]]; then
    log "taxonomy.qza existe deja : ${TAXONOMY}"
    log "Suppression de l'ancien resultat avant reclassification."
    rm -f "${TAXONOMY}"
fi

log "Classification de $(basename "${REP_SEQS}") avec SILVA 138.2"

conda run -n "${QIIME_ENV}" qiime feature-classifier classify-sklearn \
    --i-classifier "${CLASSIFIER}" \
    --i-reads "${REP_SEQS}" \
    --p-n-jobs "${NTHREADS}" \
    --o-classification "${TAXONOMY}"

[[ -f "${TAXONOMY}" ]] \
    || die "taxonomy.qza n'a pas ete produit : ${TAXONOMY}"

log "Classification reussie."
log "Resultat : ${TAXONOMY}"
