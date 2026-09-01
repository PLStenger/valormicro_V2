#!/usr/bin/env bash
# =============================================================================
# RAREFACTION + EXPORT TSV DES RESULTATS QIIME2
#
# Etape 1 : rarefaction explicite de la table ASV a SAMPLING_DEPTH reads/sample
# Etape 2 : export de la table rarefiee QZA -> BIOM -> TSV
# Etape 3 : export des metriques alpha deja calculees
#            Faith PD, Shannon, Pielou evenness, observed features
# Etape 4 : export de la taxonomie
# Etape 5 : export des matrices beta, PCoA et visualisations Emperor
#
# Attention :
# - La rarefaction est aleatoire. Ne relance pas ce script si tu veux conserver
#   exactement la meme table rarefiee, sauf si RUN_RAREFY=true est conserve.
# - Par defaut, le script ecrase la table rarefiee existante afin de la refaire.
# =============================================================================

set -Eeuo pipefail
shopt -s nullglob

# ==================== CONFIGURATION ====================

ROOTDIR="/nvme/bio/data_fungi/valormicro_V2"
RESULTS_DIR="${ROOTDIR}/02_amplicon_pipeline"

QIIME_DIR="${RESULTS_DIR}/05_qiime2"
QIIME_CORE_DIR="${QIIME_DIR}/core"
QIIME_TREE_DIR="${QIIME_DIR}/tree"
QIIME_EXPORT_DIR="${QIIME_DIR}/export"
LOG_DIR="${RESULTS_DIR}/logs"

QIIME_ENV="qiime2-amplicon-2025.7"
BIOM_ENV="biom-format"

# Profondeur de rarefaction : chaque echantillon retenu aura exactement
# SAMPLING_DEPTH reads dans la table rarefiee.
SAMPLING_DEPTH=3547

# Table source.
# Elle doit etre compatible avec tree/rooted-tree.qza, construit a partir des
# ASV presents dans au moins deux echantillons.
INPUT_TABLE="${QIIME_CORE_DIR}/table.qza"

# Arbre conserve comme controle de coherence et pour les analyses phylogenetiques.
ROOTED_TREE="${QIIME_TREE_DIR}/rooted-tree.qza"

# Artefact produit par la rarefaction explicite.
RAREFIED_TABLE="${QIIME_CORE_DIR}/RarTable_depth${SAMPLING_DEPTH}.qza"

# Artefacts alpha existants, issus du precedent core-metrics-phylogenetic.
CORE_METRICS_DIR="${QIIME_CORE_DIR}/core-metrics-depth${SAMPLING_DEPTH}"
FAITH_PD="${CORE_METRICS_DIR}/faith_pd_vector.qza"
SHANNON="${CORE_METRICS_DIR}/shannon_vector.qza"
EVENNESS="${CORE_METRICS_DIR}/evenness_vector.qza"
OBSERVED_FEATURES="${CORE_METRICS_DIR}/observed_features_vector.qza"

# Taxonomie precedemment produite par classify-sklearn.
TAXONOMY="${QIIME_CORE_DIR}/taxonomy.qza"

# Repertoires d'export.
EXPORT_CORE_DIR="${QIIME_EXPORT_DIR}/core-metrics-depth${SAMPLING_DEPTH}"
EXPORT_TAXONOMY_DIR="${QIIME_EXPORT_DIR}/taxonomy"
EXPORT_VISUAL_DIR="${QIIME_EXPORT_DIR}/visual/core-metrics-depth${SAMPLING_DEPTH}"

# true : execute qiime feature-table rarefy.
# false : reutilise RarTable_depth${SAMPLING_DEPTH}.qza deja existant.
RUN_RAREFY=true

mkdir -p \
    "${LOG_DIR}" \
    "${QIIME_CORE_DIR}" \
    "${QIIME_EXPORT_DIR}" \
    "${EXPORT_CORE_DIR}" \
    "${EXPORT_TAXONOMY_DIR}" \
    "${EXPORT_VISUAL_DIR}"

# ==================== FONCTIONS ====================

log() {
    printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOG_DIR}/04_rarefaction_export.log"
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

    [[ -f "${artifact}" ]] \
        || die "Artefact QIIME2 absent : ${artifact}"

    rm -rf "${destination}"
    mkdir -p "${destination}"

    log "Export QIIME2 : $(basename "${artifact}")"

    qiime_run tools export \
        --input-path "${artifact}" \
        --output-path "${destination}"
}

# ==================== CONTROLES ====================

log "Debut de la rarefaction et des exports QIIME2"
log "Environnement QIIME2 : ${QIIME_ENV}"
log "Environnement BIOM : ${BIOM_ENV}"
log "Profondeur de rarefaction : ${SAMPLING_DEPTH}"

command -v conda >/dev/null 2>&1 \
    || die "Conda est introuvable dans le PATH."

qiime_run --version \
    || die "QIIME2 ne demarre pas dans l'environnement ${QIIME_ENV}."

conda run -n "${BIOM_ENV}" biom --help >/dev/null \
    || die "La commande biom est absente de l'environnement ${BIOM_ENV}."

[[ -f "${INPUT_TABLE}" ]] \
    || die "Table ASV source absente : ${INPUT_TABLE}"

[[ -f "${ROOTED_TREE}" ]] \
    || die "Arbre enracine absent : ${ROOTED_TREE}"

[[ -f "${TAXONOMY}" ]] \
    || die "Taxonomie absente : ${TAXONOMY}"

# ==================== 1. RAREFACTION EXPLICITE ====================

if [[ "${RUN_RAREFY}" == true ]]; then
    if [[ -f "${RAREFIED_TABLE}" ]]; then
        log "Suppression de l'ancienne table rarefiee : ${RAREFIED_TABLE}"
        rm -f "${RAREFIED_TABLE}"
    fi

    log "Rarefaction de la table ASV a ${SAMPLING_DEPTH} reads par echantillon"
    log "Table source : ${INPUT_TABLE}"

    qiime_run feature-table rarefy \
        --i-table "${INPUT_TABLE}" \
        --p-sampling-depth "${SAMPLING_DEPTH}" \
        --p-no-with-replacement \
        --o-rarefied-table "${RAREFIED_TABLE}"
else
    log "Rarefaction desactivee : reutilisation de ${RAREFIED_TABLE}"
fi

[[ -f "${RAREFIED_TABLE}" ]] \
    || die "Table rarefiee absente apres rarefaction : ${RAREFIED_TABLE}"

log "Table rarefiee QIIME2 disponible : ${RAREFIED_TABLE}"

# ==================== 2. EXPORT TABLE RAREFIEE ====================

TABLE_EXPORT_DIR="${EXPORT_CORE_DIR}/rarefied_table"
BIOM_TABLE="${TABLE_EXPORT_DIR}/feature-table.biom"
TSV_TABLE="${TABLE_EXPORT_DIR}/table-from-biom.tsv"
ASV_TABLE="${TABLE_EXPORT_DIR}/ASV.tsv"

export_qza "${RAREFIED_TABLE}" "${TABLE_EXPORT_DIR}"

[[ -f "${BIOM_TABLE}" ]] \
    || die "Fichier BIOM absent apres export : ${BIOM_TABLE}"

log "Conversion de la table BIOM vers TSV"

conda run -n "${BIOM_ENV}" biom convert \
    -i "${BIOM_TABLE}" \
    -o "${TSV_TABLE}" \
    --to-tsv

# biom convert produit :
# ligne 1 : commentaire de construction BIOM
# ligne 2 : #OTU ID <tab> sample_1 <tab> ...
# On retire uniquement la ligne commentaire et renomme #OTU ID en ASV_ID.
sed '1d; s/^#OTU ID/ASV_ID/' "${TSV_TABLE}" > "${ASV_TABLE}"

[[ -s "${ASV_TABLE}" ]] \
    || die "Echec de creation de la table ASV TSV : ${ASV_TABLE}"

log "Table ASV rarefiee TSV produite : ${ASV_TABLE}"

# ==================== 3. METRIQUES ALPHA ====================

# Ces quatre artefacts doivent deja avoir ete produits par :
# qiime diversity core-metrics-phylogenetic
#
# Ils sont seulement exportes ici. Ce script ne recalcule pas Faith PD,
# Shannon, evenness ou observed features.

[[ -f "${FAITH_PD}" ]] \
    || die "Faith PD absent : ${FAITH_PD}"

[[ -f "${SHANNON}" ]] \
    || die "Shannon absent : ${SHANNON}"

[[ -f "${EVENNESS}" ]] \
    || die "Evenness absent : ${EVENNESS}"

[[ -f "${OBSERVED_FEATURES}" ]] \
    || die "Observed features absent : ${OBSERVED_FEATURES}"

log "Export des metriques alpha-diversite"

export_qza "${FAITH_PD}" "${EXPORT_CORE_DIR}/faith_pd"
export_qza "${SHANNON}" "${EXPORT_CORE_DIR}/shannon"
export_qza "${EVENNESS}" "${EXPORT_CORE_DIR}/evenness"
export_qza "${OBSERVED_FEATURES}" "${EXPORT_CORE_DIR}/observed_features"

# ==================== 4. TAXONOMIE ====================

log "Export de la taxonomie"

export_qza "${TAXONOMY}" "${EXPORT_TAXONOMY_DIR}"

[[ -f "${EXPORT_TAXONOMY_DIR}/taxonomy.tsv" ]] \
    || die "taxonomy.tsv absent apres export : ${EXPORT_TAXONOMY_DIR}/taxonomy.tsv"

# ==================== 5. BETA DIVERSITE ====================

# Les resultats beta doivent deja exister depuis core-metrics-phylogenetic.
# Ils sont exportes ici au format texte utilisable hors QIIME2.

log "Export des distances beta et des resultats PCoA"

for artifact in \
    "${CORE_METRICS_DIR}/bray_curtis_distance_matrix.qza" \
    "${CORE_METRICS_DIR}/jaccard_distance_matrix.qza" \
    "${CORE_METRICS_DIR}/unweighted_unifrac_distance_matrix.qza" \
    "${CORE_METRICS_DIR}/weighted_unifrac_distance_matrix.qza" \
    "${CORE_METRICS_DIR}/bray_curtis_pcoa_results.qza" \
    "${CORE_METRICS_DIR}/jaccard_pcoa_results.qza" \
    "${CORE_METRICS_DIR}/unweighted_unifrac_pcoa_results.qza" \
    "${CORE_METRICS_DIR}/weighted_unifrac_pcoa_results.qza"; do

    [[ -f "${artifact}" ]] \
        || die "Artefact beta-diversite absent : ${artifact}"

    artifact_name="$(basename "${artifact}" .qza)"

    export_qza \
        "${artifact}" \
        "${EXPORT_CORE_DIR}/${artifact_name}"
done

# ==================== 6. EMPEROR QZV ====================

log "Export des visualisations Emperor"

for qzv in \
    "${CORE_METRICS_DIR}/bray_curtis_emperor.qzv" \
    "${CORE_METRICS_DIR}/jaccard_emperor.qzv" \
    "${CORE_METRICS_DIR}/unweighted_unifrac_emperor.qzv" \
    "${CORE_METRICS_DIR}/weighted_unifrac_emperor.qzv"; do

    [[ -f "${qzv}" ]] \
        || die "Visualisation Emperor absente : ${qzv}"

    qzv_name="$(basename "${qzv}" .qzv)"

    export_qza \
        "${qzv}" \
        "${EXPORT_VISUAL_DIR}/${qzv_name}"
done

# ==================== FIN ====================

log "Pipeline termine avec succes."
log "Table QIIME2 rarefiee : ${RAREFIED_TABLE}"
log "Table ASV TSV rarefiee : ${ASV_TABLE}"
log "Faith PD TSV : ${EXPORT_CORE_DIR}/faith_pd/alpha-diversity.tsv"
log "Shannon TSV : ${EXPORT_CORE_DIR}/shannon/alpha-diversity.tsv"
log "Evenness TSV : ${EXPORT_CORE_DIR}/evenness/alpha-diversity.tsv"
log "Observed features TSV : ${EXPORT_CORE_DIR}/observed_features/alpha-diversity.tsv"
log "Taxonomie TSV : ${EXPORT_TAXONOMY_DIR}/taxonomy.tsv"
log "Exports beta-diversite : ${EXPORT_CORE_DIR}"
log "Exports Emperor : ${EXPORT_VISUAL_DIR}"
