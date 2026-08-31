#!/usr/bin/env bash
# =============================================================================
# Reprise : export TSV des resultats QIIME2 deja calcules
# - table ASV rarefiee : BIOM -> TSV
# - metriques alpha : Faith PD, Shannon, evenness, observed features
# - taxonomie
# - distances beta et PCoA
# - visualisations QZV
# =============================================================================

set -Eeuo pipefail
shopt -s nullglob

# ==================== CONFIGURATION ====================

ROOTDIR="/nvme/bio/data_fungi/valormicro_V2"
RESULTS_DIR="${ROOTDIR}/02_amplicon_pipeline"

QIIME_DIR="${RESULTS_DIR}/05_qiime2"
QIIME_CORE_DIR="${QIIME_DIR}/core"
QIIME_EXPORT_DIR="${QIIME_DIR}/export"
LOG_DIR="${RESULTS_DIR}/logs"

QIIME_ENV="qiime2-amplicon-2025.7"
BIOM_ENV="biom-format"

SAMPLING_DEPTH=3547
CORE_METRICS_DIR="${QIIME_CORE_DIR}/core-metrics-depth${SAMPLING_DEPTH}"

RAREFIED_TABLE="${CORE_METRICS_DIR}/rarefied_table.qza"
FAITH_PD="${CORE_METRICS_DIR}/faith_pd_vector.qza"
SHANNON="${CORE_METRICS_DIR}/shannon_vector.qza"
EVENNESS="${CORE_METRICS_DIR}/evenness_vector.qza"
OBSERVED_FEATURES="${CORE_METRICS_DIR}/observed_features_vector.qza"

TAXONOMY="${QIIME_CORE_DIR}/taxonomy.qza"

EXPORT_CORE_DIR="${QIIME_EXPORT_DIR}/core-metrics-depth${SAMPLING_DEPTH}"
EXPORT_TAXONOMY_DIR="${QIIME_EXPORT_DIR}/taxonomy"
EXPORT_VISUAL_DIR="${QIIME_EXPORT_DIR}/visual/core-metrics-depth${SAMPLING_DEPTH}"

mkdir -p "${LOG_DIR}" "${EXPORT_CORE_DIR}" "${EXPORT_TAXONOMY_DIR}" "${EXPORT_VISUAL_DIR}"

# ==================== FONCTIONS ====================

log() {
    printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOG_DIR}/04_export_core_metrics.log"
}

die() {
    log "ERREUR : $*"
    exit 1
}

trap 'rc=$?; log "ERREUR : code ${rc}, ligne ${LINENO} : ${BASH_COMMAND}"; exit "${rc}"' ERR

qiime_run() {
    conda run -n "${QIIME_ENV}" qiime "$@"
}

export_qza() {
    local artifact="$1"
    local destination="$2"

    [[ -f "${artifact}" ]] || die "Artefact QIIME2 absent : ${artifact}"

    rm -rf "${destination}"
    mkdir -p "${destination}"

    log "Export : $(basename "${artifact}")"
    qiime_run tools export \
        --input-path "${artifact}" \
        --output-path "${destination}"
}

# ==================== CONTROLES ====================

log "Debut des exports des resultats deja calcules"

command -v conda >/dev/null 2>&1 \
    || die "Conda est introuvable dans le PATH."

qiime_run --version \
    || die "QIIME2 ne demarre pas dans ${QIIME_ENV}."

conda run -n "${BIOM_ENV}" biom --help >/dev/null \
    || die "La commande biom est absente de l'environnement ${BIOM_ENV}."

[[ -f "${RAREFIED_TABLE}" ]] \
    || die "Table rarefiee absente : ${RAREFIED_TABLE}"

[[ -f "${FAITH_PD}" ]] \
    || die "Faith PD absent : ${FAITH_PD}"

[[ -f "${SHANNON}" ]] \
    || die "Shannon absent : ${SHANNON}"

[[ -f "${EVENNESS}" ]] \
    || die "Evenness absent : ${EVENNESS}"

[[ -f "${OBSERVED_FEATURES}" ]] \
    || die "Observed features absent : ${OBSERVED_FEATURES}"

[[ -f "${TAXONOMY}" ]] \
    || die "Taxonomie absente : ${TAXONOMY}"

# ==================== TABLE ASV RAREFIEE ====================

TABLE_EXPORT_DIR="${EXPORT_CORE_DIR}/rarefied_table"

export_qza "${RAREFIED_TABLE}" "${TABLE_EXPORT_DIR}"

BIOM_TABLE="${TABLE_EXPORT_DIR}/feature-table.biom"
TSV_TABLE="${TABLE_EXPORT_DIR}/table-from-biom.tsv"
ASV_TABLE="${TABLE_EXPORT_DIR}/ASV.tsv"

[[ -f "${BIOM_TABLE}" ]] \
    || die "Fichier BIOM absent apres export : ${BIOM_TABLE}"

log "Conversion BIOM vers TSV"

conda run -n "${BIOM_ENV}" biom convert \
    -i "${BIOM_TABLE}" \
    -o "${TSV_TABLE}" \
    --to-tsv

# biom convert ajoute une ligne de commentaire ; elle est retiree ici.
# L'en-tete #OTU ID est renomme ASV_ID.
sed '1d; s/^#OTU ID/ASV_ID/' "${TSV_TABLE}" > "${ASV_TABLE}"

[[ -s "${ASV_TABLE}" ]] \
    || die "Echec de creation de la table ASV TSV : ${ASV_TABLE}"

# ==================== METRIQUES ALPHA ====================

export_qza "${FAITH_PD}" "${EXPORT_CORE_DIR}/faith_pd"
export_qza "${SHANNON}" "${EXPORT_CORE_DIR}/shannon"
export_qza "${EVENNESS}" "${EXPORT_CORE_DIR}/evenness"
export_qza "${OBSERVED_FEATURES}" "${EXPORT_CORE_DIR}/observed_features"

# ==================== TAXONOMIE ====================

export_qza "${TAXONOMY}" "${EXPORT_TAXONOMY_DIR}"

[[ -f "${EXPORT_TAXONOMY_DIR}/taxonomy.tsv" ]] \
    || die "taxonomy.tsv absent apres export : ${EXPORT_TAXONOMY_DIR}/taxonomy.tsv"

# ==================== BETA DIVERSITE ====================

for artifact in \
    "${CORE_METRICS_DIR}/bray_curtis_distance_matrix.qza" \
    "${CORE_METRICS_DIR}/jaccard_distance_matrix.qza" \
    "${CORE_METRICS_DIR}/unweighted_unifrac_distance_matrix.qza" \
    "${CORE_METRICS_DIR}/weighted_unifrac_distance_matrix.qza" \
    "${CORE_METRICS_DIR}/bray_curtis_pcoa_results.qza" \
    "${CORE_METRICS_DIR}/jaccard_pcoa_results.qza" \
    "${CORE_METRICS_DIR}/unweighted_unifrac_pcoa_results.qza" \
    "${CORE_METRICS_DIR}/weighted_unifrac_pcoa_results.qza"; do

    artifact_name="$(basename "${artifact}" .qza)"
    export_qza "${artifact}" "${EXPORT_CORE_DIR}/${artifact_name}"
done

# ==================== VISUALISATIONS ====================

for qzv in \
    "${CORE_METRICS_DIR}/bray_curtis_emperor.qzv" \
    "${CORE_METRICS_DIR}/jaccard_emperor.qzv" \
    "${CORE_METRICS_DIR}/unweighted_unifrac_emperor.qzv" \
    "${CORE_METRICS_DIR}/weighted_unifrac_emperor.qzv"; do

    [[ -f "${qzv}" ]] || die "Visualisation absente : ${qzv}"

    qzv_name="$(basename "${qzv}" .qzv)"
    export_qza "${qzv}" "${EXPORT_VISUAL_DIR}/${qzv_name}"
done

# ==================== FIN ====================

log "Exports termines avec succes."
log "Table ASV TSV : ${ASV_TABLE}"
log "Faith PD TSV : ${EXPORT_CORE_DIR}/faith_pd/alpha-diversity.tsv"
log "Shannon TSV : ${EXPORT_CORE_DIR}/shannon/alpha-diversity.tsv"
log "Evenness TSV : ${EXPORT_CORE_DIR}/evenness/alpha-diversity.tsv"
log "Observed features TSV : ${EXPORT_CORE_DIR}/observed_features/alpha-diversity.tsv"
log "Taxonomie TSV : ${EXPORT_TAXONOMY_DIR}/taxonomy.tsv"
