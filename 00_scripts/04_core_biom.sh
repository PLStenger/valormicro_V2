#!/usr/bin/env bash
# =============================================================================
# Core diversity metrics + core microbiome
# QIIME2 amplicon 2025.7
# =============================================================================

set -Eeuo pipefail
shopt -s nullglob

# ==================== CONFIGURATION ====================

ROOTDIR="/nvme/bio/data_fungi/valormicro_V2"
PROJECT_DIR="${ROOTDIR}"
RESULTS_DIR="${PROJECT_DIR}/02_amplicon_pipeline"

DATABASE="${RESULTS_DIR}/04_database_files"
QIIME_DIR="${RESULTS_DIR}/05_qiime2"
QIIME_CORE_DIR="${QIIME_DIR}/core"
QIIME_TREE_DIR="${QIIME_DIR}/tree"
QIIME_VISUAL_DIR="${QIIME_DIR}/visual"
QIIME_EXPORT_DIR="${QIIME_DIR}/export"
LOG_DIR="${RESULTS_DIR}/logs"

NTHREADS=16
SAMPLING_DEPTH=3547

TMPDIR="${ROOTDIR}/tmp"
QIIME_ENV="qiime2-amplicon-2025.7"

# La table est celle compatible avec l'arbre, construit depuis les ASV
# conservés dans au moins deux échantillons.
TABLE="${QIIME_CORE_DIR}/table_min2samples.qza"
ROOTED_TREE="${QIIME_TREE_DIR}/rooted-tree.qza"
METADATA="${DATABASE}/sample-metadata.tsv"
TAXONOMY="${QIIME_CORE_DIR}/taxonomy.qza"

# Toutes les sorties core-metrics sont centralisées ici.
CORE_METRICS_DIR="${QIIME_CORE_DIR}/core-metrics-depth${SAMPLING_DEPTH}"

# Exports lisibles hors artefacts QIIME2.
EXPORT_CORE_DIR="${QIIME_EXPORT_DIR}/core-metrics-depth${SAMPLING_DEPTH}"
EXPORT_VISUAL_DIR="${QIIME_EXPORT_DIR}/visual/core-metrics-depth${SAMPLING_DEPTH}"
EXPORT_TAXONOMY_DIR="${QIIME_EXPORT_DIR}/taxonomy"

mkdir -p \
    "${TMPDIR}" \
    "${LOG_DIR}" \
    "${QIIME_EXPORT_DIR}" \
    "${EXPORT_CORE_DIR}" \
    "${EXPORT_VISUAL_DIR}" \
    "${EXPORT_TAXONOMY_DIR}"

export TMPDIR

# ==================== FONCTIONS ====================

log() {
    printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOG_DIR}/04_core_biom.log"
}

die() {
    log "ERREUR : $*"
    exit 1
}

trap 'rc=$?; log "ERREUR : code ${rc}, ligne ${LINENO} : ${BASH_COMMAND}"; exit "${rc}"' ERR

qiime_run() {
    conda run -n "${QIIME_ENV}" qiime "$@"
}

# ==================== CONTROLES ====================

log "Debut du calcul des core diversity metrics"
log "Environnement QIIME2 : ${QIIME_ENV}"
log "Nombre de threads : ${NTHREADS}"
log "Profondeur de rarefaction : ${SAMPLING_DEPTH}"

command -v conda >/dev/null 2>&1 \
    || die "Conda est introuvable dans le PATH."

qiime_run --version \
    || die "QIIME2 ne demarre pas dans l'environnement ${QIIME_ENV}."

[[ -f "${TABLE}" ]] \
    || die "Table ASV compatible avec l'arbre absente : ${TABLE}"

[[ -f "${ROOTED_TREE}" ]] \
    || die "Arbre enracine absent : ${ROOTED_TREE}"

[[ -f "${METADATA}" ]] \
    || die "Metadata absente : ${METADATA}"

[[ -f "${TAXONOMY}" ]] \
    || die "Taxonomie absente : ${TAXONOMY}"

# ==================== CORE METRICS ====================

cd "${QIIME_DIR}" || die "Impossible d'acceder a ${QIIME_DIR}"

if [[ -d "${CORE_METRICS_DIR}" ]]; then
    log "Suppression de l'ancien repertoire core metrics : ${CORE_METRICS_DIR}"
    rm -rf "${CORE_METRICS_DIR}"
fi

log "Calcul des metriques alpha et beta de diversite"

qiime_run diversity core-metrics-phylogenetic \
    --i-phylogeny "${ROOTED_TREE}" \
    --i-table "${TABLE}" \
    --p-sampling-depth "${SAMPLING_DEPTH}" \
    --m-metadata-file "${METADATA}" \
    --p-n-jobs-or-threads "${NTHREADS}" \
    --output-dir "${CORE_METRICS_DIR}"

RAREFIED_TABLE="${CORE_METRICS_DIR}/rarefied_table.qza"
FAITH_PD="${CORE_METRICS_DIR}/faith_pd_vector.qza"
EVENNESS="${CORE_METRICS_DIR}/evenness_vector.qza"
SHANNON="${CORE_METRICS_DIR}/shannon_vector.qza"
OBSERVED_ASV="${CORE_METRICS_DIR}/observed_features_vector.qza"

[[ -f "${RAREFIED_TABLE}" ]] \
    || die "Table raréfiee absente : ${RAREFIED_TABLE}"

[[ -f "${FAITH_PD}" ]] \
    || die "Vecteur Faith PD absent : ${FAITH_PD}"

[[ -f "${EVENNESS}" ]] \
    || die "Vecteur Pielou absent : ${EVENNESS}"

[[ -f "${SHANNON}" ]] \
    || die "Vecteur Shannon absent : ${SHANNON}"

[[ -f "${OBSERVED_ASV}" ]] \
    || die "Vecteur observed features absent : ${OBSERVED_ASV}"

log "Core diversity metrics terminees"

# ==================== CORE MICROBIOME ====================

CORE_BIOM_QZV="${QIIME_VISUAL_DIR}/CoreBiom-all-depth${SAMPLING_DEPTH}.qzv"

rm -f "${CORE_BIOM_QZV}"

log "Calcul des ASV du core microbiome"

qiime_run feature-table core-features \
    --i-table "${RAREFIED_TABLE}" \
    --p-min-fraction 0.1 \
    --p-max-fraction 1.0 \
    --p-steps 10 \
    --o-visualization "${CORE_BIOM_QZV}"

# ==================== EXPORT TABLE RAREFIEE ====================

rm -rf "${EXPORT_CORE_DIR}/rarefied_table"
mkdir -p "${EXPORT_CORE_DIR}/rarefied_table"

log "Export de la table rarefiee"

qiime_run tools export \
    --input-path "${RAREFIED_TABLE}" \
    --output-path "${EXPORT_CORE_DIR}/rarefied_table"

BIOM_TABLE="${EXPORT_CORE_DIR}/rarefied_table/feature-table.biom"
TSV_TABLE="${EXPORT_CORE_DIR}/rarefied_table/table-from-biom.tsv"
ASV_TABLE="${EXPORT_CORE_DIR}/rarefied_table/ASV.tsv"

[[ -f "${BIOM_TABLE}" ]] \
    || die "Export BIOM absent : ${BIOM_TABLE}"

command -v biom >/dev/null 2>&1 \
    || die "La commande biom est introuvable. Activez ou installez l'environnement biom-format."

biom convert \
    -i "${BIOM_TABLE}" \
    -o "${TSV_TABLE}" \
    --to-tsv

sed '1d; s/^#OTU ID/ASV_ID/' "${TSV_TABLE}" > "${ASV_TABLE}"

# ==================== EXPORT METRIQUES ALPHA ====================

log "Export des vecteurs alpha-diversite"

for artifact in \
    "${FAITH_PD}" \
    "${EVENNESS}" \
    "${SHANNON}" \
    "${OBSERVED_ASV}"; do

    artifact_name="$(basename "${artifact}" .qza)"
    artifact_export_dir="${EXPORT_CORE_DIR}/${artifact_name}"

    rm -rf "${artifact_export_dir}"
    mkdir -p "${artifact_export_dir}"

    qiime_run tools export \
        --input-path "${artifact}" \
        --output-path "${artifact_export_dir}"
done

# ==================== EXPORT METRIQUES BETA ====================

log "Export des matrices de distance et PCoA"

for artifact in \
    "${CORE_METRICS_DIR}/unweighted_unifrac_distance_matrix.qza" \
    "${CORE_METRICS_DIR}/weighted_unifrac_distance_matrix.qza" \
    "${CORE_METRICS_DIR}/jaccard_distance_matrix.qza" \
    "${CORE_METRICS_DIR}/bray_curtis_distance_matrix.qza" \
    "${CORE_METRICS_DIR}/unweighted_unifrac_pcoa_results.qza" \
    "${CORE_METRICS_DIR}/weighted_unifrac_pcoa_results.qza" \
    "${CORE_METRICS_DIR}/jaccard_pcoa_results.qza" \
    "${CORE_METRICS_DIR}/bray_curtis_pcoa_results.qza"; do

    [[ -f "${artifact}" ]] || die "Artefact beta-diversite absent : ${artifact}"

    artifact_name="$(basename "${artifact}" .qza)"
    artifact_export_dir="${EXPORT_CORE_DIR}/${artifact_name}"

    rm -rf "${artifact_export_dir}"
    mkdir -p "${artifact_export_dir}"

    qiime_run tools export \
        --input-path "${artifact}" \
        --output-path "${artifact_export_dir}"
done

# ==================== EXPORT VISUALISATIONS ====================

log "Export des visualisations QIIME2"

for qzv in \
    "${CORE_BIOM_QZV}" \
    "${CORE_METRICS_DIR}/unweighted_unifrac_emperor.qzv" \
    "${CORE_METRICS_DIR}/weighted_unifrac_emperor.qzv" \
    "${CORE_METRICS_DIR}/jaccard_emperor.qzv" \
    "${CORE_METRICS_DIR}/bray_curtis_emperor.qzv"; do

    [[ -f "${qzv}" ]] || die "Visualisation absente : ${qzv}"

    qzv_name="$(basename "${qzv}" .qzv)"
    qzv_export_dir="${EXPORT_VISUAL_DIR}/${qzv_name}"

    rm -rf "${qzv_export_dir}"
    mkdir -p "${qzv_export_dir}"

    qiime_run tools export \
        --input-path "${qzv}" \
        --output-path "${qzv_export_dir}"
done

# ==================== EXPORT TAXONOMIE ====================

rm -rf "${EXPORT_TAXONOMY_DIR}"
mkdir -p "${EXPORT_TAXONOMY_DIR}"

log "Export de la taxonomie"

qiime_run tools export \
    --input-path "${TAXONOMY}" \
    --output-path "${EXPORT_TAXONOMY_DIR}"

[[ -f "${EXPORT_TAXONOMY_DIR}/taxonomy.tsv" ]] \
    || die "taxonomy.tsv absent apres export."

log "Pipeline core biom termine avec succes."
log "Resultats QIIME2 : ${CORE_METRICS_DIR}"
log "Exports lisibles : ${EXPORT_CORE_DIR}"
log "Taxonomie TSV : ${EXPORT_TAXONOMY_DIR}/taxonomy.tsv"
