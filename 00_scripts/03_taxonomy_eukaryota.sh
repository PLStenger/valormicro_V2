#!/usr/bin/env bash
# =============================================================================
# Classification des ASV eucaryotes avec PR2 v5
#
# Projet : valormicro_V2
# Strategie :
# 1. Classifier tous les ASV avec SILVA 138.2 (resultat deja produit par
#    03_taxonomy_16S.sh, reutilise s'il est present).
# 2. Retirer les ASV annotés Bacteria, Archaea, mitochondries et chloroplastes
#    par SILVA afin d'obtenir un jeu de candidats eucaryotes.
# 3. Classifier ces candidats avec le classifieur PR2 v5 515F-926R existant.
# 4. Produire les artefacts QZA/QZV et les exports TSV/FASTA/BIOM.
#
# Le script n'ecrase pas les fichiers existants sauf si FORCE=true.
# =============================================================================

set -Eeuo pipefail
shopt -s nullglob

# ==================== CONFIGURATION ====================

ROOTDIR="/nvme/bio/data_fungi/valormicro_V2"
RESULTS_DIR="${ROOTDIR}/02_amplicon_pipeline"

QIIME_ENV="qiime2-amplicon-2025.7"
BIOM_ENV="biom-format"
NTHREADS=16
TMPDIR="${ROOTDIR}/tmp"

QIIME_DIR="${RESULTS_DIR}/05_qiime2"
QIIME_CORE_DIR="${QIIME_DIR}/core"
QIIME_VISUAL_DIR="${QIIME_DIR}/visual"
QIIME_EXPORT_DIR="${QIIME_DIR}/export"
DATABASE_DIR="${RESULTS_DIR}/04_database_files"
LOG_DIR="${RESULTS_DIR}/logs"

# Entrees du pipeline amplicon deja execute.
TABLE="${QIIME_CORE_DIR}/table.qza"
REP_SEQS="${QIIME_CORE_DIR}/rep-seqs.qza"
SILVA_TAXONOMY="${QIIME_CORE_DIR}/taxonomy.qza"
METADATA="${DATABASE_DIR}/sample-metadata.tsv"

# Classifieur PR2 deja construit dans l'ancien projet.
# Ce classifieur est base sur PR2 v5, avec extraction in silico 515F-926R.
PR2_CLASSIFIER="/nvme/bio/data_fungi/vague_project/98_databasefiles/pr2-v5-classifier-515F-926R.qza"

# Sorties eucaryotes dediees : aucun fichier 16S/SILVA n'est ecrase.
EUK_DIR="${QIIME_CORE_DIR}/eukaryota"
EUK_TABLE_TMP="${EUK_DIR}/table-eukaryota-no-prokaryotes.qza"
EUK_TABLE="${EUK_DIR}/table-eukaryota.qza"
EUK_REP_SEQS_TMP="${EUK_DIR}/rep-seqs-eukaryota-no-prokaryotes.qza"
EUK_REP_SEQS="${EUK_DIR}/rep-seqs-eukaryota.qza"
EUK_TAXONOMY="${EUK_DIR}/taxonomy-pr2-eukaryota.qza"

EUK_TABLE_QZV="${QIIME_VISUAL_DIR}/table-eukaryota-summary.qzv"
EUK_TAXONOMY_QZV="${QIIME_VISUAL_DIR}/taxonomy-pr2-eukaryota.qzv"
EUK_BAR_PLOT_QZV="${QIIME_VISUAL_DIR}/taxa-bar-plots-eukaryota-pr2.qzv"

EUK_EXPORT_DIR="${QIIME_EXPORT_DIR}/eukaryota"
EUK_TABLE_EXPORT_DIR="${EUK_EXPORT_DIR}/table"
EUK_REP_SEQS_EXPORT_DIR="${EUK_EXPORT_DIR}/rep-seqs"
EUK_TAXONOMY_EXPORT_DIR="${EUK_EXPORT_DIR}/taxonomy-pr2"
EUK_ASV_TSV="${EUK_TABLE_EXPORT_DIR}/ASV.tsv"
EUK_ASV_TAXONOMY_TSV="${EUK_EXPORT_DIR}/ASV_taxonomy_PR2.tsv"

# true : recalcule et ecrase les sorties eucaryotes existantes.
# false : reutilise chaque sortie deja produite lorsqu'elle existe.
FORCE=false

mkdir -p \
    "${TMPDIR}" \
    "${LOG_DIR}" \
    "${EUK_DIR}" \
    "${QIIME_VISUAL_DIR}" \
    "${EUK_EXPORT_DIR}"

export TMPDIR

# ==================== FONCTIONS ====================

log() {
    printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOG_DIR}/03_taxonomy_eukaryota.log"
}

die() {
    log "ERREUR : $*"
    exit 1
}

trap 'rc=$?; log "ERREUR : code ${rc}, ligne ${LINENO} : ${BASH_COMMAND}"; exit "${rc}"' ERR

qiime_run() {
    conda run -n "${QIIME_ENV}" qiime "$@"
}

remove_output_if_forced() {
    local output="$1"
    if [[ "${FORCE}" == true && -e "${output}" ]]; then
        log "Suppression de la sortie existante : ${output}"
        rm -rf "${output}"
    fi
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

log "Debut de la classification des eucaryotes avec PR2 v5"
log "Environnement QIIME2 : ${QIIME_ENV}"
log "Nombre de jobs classify-sklearn : ${NTHREADS}"
log "FORCE : ${FORCE}"

command -v conda >/dev/null 2>&1 || die "Conda est introuvable dans le PATH."

qiime_run --version || die "QIIME2 ne demarre pas dans l'environnement ${QIIME_ENV}."

for required_file in \
    "${TABLE}" \
    "${REP_SEQS}" \
    "${SILVA_TAXONOMY}" \
    "${METADATA}" \
    "${PR2_CLASSIFIER}"; do
    [[ -f "${required_file}" ]] || die "Fichier requis absent : ${required_file}"
done

log "Validation du classifieur PR2"
qiime_run tools validate "${PR2_CLASSIFIER}" \
    || die "Classifieur PR2 invalide : ${PR2_CLASSIFIER}"

# ==================== CANDIDATS EUCARYOTES ====================

# SILVA classe ici tous les ASV, y compris les reads non bacteriens captures par
# les amorces 515F-926R. On retire d'abord les taxons prokaryotes sans utiliser
# de filtre include : les ASV restant forment le jeu de candidats eucaryotes.

remove_output_if_forced "${EUK_TABLE_TMP}"
if [[ ! -f "${EUK_TABLE_TMP}" ]]; then
    log "Filtrage de Bacteria et Archaea dans la table"
    qiime_run taxa filter-table \
        --i-table "${TABLE}" \
        --i-taxonomy "${SILVA_TAXONOMY}" \
        --p-mode contains \
        --p-exclude "D_0__Bacteria,D_0__Archaea,k__Bacteria,k__Archaea" \
        --o-filtered-table "${EUK_TABLE_TMP}"
else
    log "Table intermediaire eucaryote deja presente : ${EUK_TABLE_TMP}"
fi

remove_output_if_forced "${EUK_TABLE}"
if [[ ! -f "${EUK_TABLE}" ]]; then
    log "Retrait des mitochondries et chloroplastes de la table"
    qiime_run taxa filter-table \
        --i-table "${EUK_TABLE_TMP}" \
        --i-taxonomy "${SILVA_TAXONOMY}" \
        --p-mode contains \
        --p-exclude "Mitochondria,Chloroplast" \
        --o-filtered-table "${EUK_TABLE}"
else
    log "Table eucaryote deja presente : ${EUK_TABLE}"
fi

remove_output_if_forced "${EUK_REP_SEQS_TMP}"
if [[ ! -f "${EUK_REP_SEQS_TMP}" ]]; then
    log "Filtrage de Bacteria et Archaea dans les sequences representatives"
    qiime_run taxa filter-seqs \
        --i-sequences "${REP_SEQS}" \
        --i-taxonomy "${SILVA_TAXONOMY}" \
        --p-mode contains \
        --p-exclude "D_0__Bacteria,D_0__Archaea,k__Bacteria,k__Archaea" \
        --o-filtered-sequences "${EUK_REP_SEQS_TMP}"
else
    log "Sequences intermediaires eucaryotes deja presentes : ${EUK_REP_SEQS_TMP}"
fi

remove_output_if_forced "${EUK_REP_SEQS}"
if [[ ! -f "${EUK_REP_SEQS}" ]]; then
    log "Retrait des mitochondries et chloroplastes des sequences representatives"
    qiime_run taxa filter-seqs \
        --i-sequences "${EUK_REP_SEQS_TMP}" \
        --i-taxonomy "${SILVA_TAXONOMY}" \
        --p-mode contains \
        --p-exclude "Mitochondria,Chloroplast" \
        --o-filtered-sequences "${EUK_REP_SEQS}"
else
    log "Sequences eucaryotes deja presentes : ${EUK_REP_SEQS}"
fi

remove_output_if_forced "${EUK_TABLE_QZV}"
if [[ ! -f "${EUK_TABLE_QZV}" ]]; then
    log "Resume de la table des candidats eucaryotes"
    qiime_run feature-table summarize \
        --i-table "${EUK_TABLE}" \
        --m-sample-metadata-file "${METADATA}" \
        --o-visualization "${EUK_TABLE_QZV}"
fi

# ==================== CLASSIFICATION PR2 ====================

remove_output_if_forced "${EUK_TAXONOMY}"
if [[ ! -f "${EUK_TAXONOMY}" ]]; then
    log "Classification PR2 v5 des candidats eucaryotes"
    qiime_run feature-classifier classify-sklearn \
        --i-classifier "${PR2_CLASSIFIER}" \
        --i-reads "${EUK_REP_SEQS}" \
        --p-n-jobs "${NTHREADS}" \
        --p-confidence 0.7 \
        --o-classification "${EUK_TAXONOMY}"
else
    log "Taxonomie PR2 deja presente : ${EUK_TAXONOMY}"
fi

[[ -f "${EUK_TAXONOMY}" ]] \
    || die "La classification PR2 n'a pas produit : ${EUK_TAXONOMY}"

remove_output_if_forced "${EUK_TAXONOMY_QZV}"
if [[ ! -f "${EUK_TAXONOMY_QZV}" ]]; then
    log "Creation de la visualisation de taxonomie PR2"
    qiime_run metadata tabulate \
        --m-input-file "${EUK_TAXONOMY}" \
        --o-visualization "${EUK_TAXONOMY_QZV}"
fi

remove_output_if_forced "${EUK_BAR_PLOT_QZV}"
if [[ ! -f "${EUK_BAR_PLOT_QZV}" ]]; then
    log "Creation du barplot taxonomique eucaryote PR2"
    qiime_run taxa barplot \
        --i-table "${EUK_TABLE}" \
        --i-taxonomy "${EUK_TAXONOMY}" \
        --m-metadata-file "${METADATA}" \
        --o-visualization "${EUK_BAR_PLOT_QZV}"
fi

# ==================== EXPORTS ====================

log "Export de la table eucaryote, des sequences et de la taxonomie PR2"

export_qza "${EUK_TABLE}" "${EUK_TABLE_EXPORT_DIR}"
export_qza "${EUK_REP_SEQS}" "${EUK_REP_SEQS_EXPORT_DIR}"
export_qza "${EUK_TAXONOMY}" "${EUK_TAXONOMY_EXPORT_DIR}"

BIOM_TABLE="${EUK_TABLE_EXPORT_DIR}/feature-table.biom"
TSV_TABLE="${EUK_TABLE_EXPORT_DIR}/table-from-biom.tsv"
TAXONOMY_TSV="${EUK_TAXONOMY_EXPORT_DIR}/taxonomy.tsv"

[[ -f "${BIOM_TABLE}" ]] || die "Export BIOM absent : ${BIOM_TABLE}"
[[ -f "${TAXONOMY_TSV}" ]] || die "Export taxonomie absent : ${TAXONOMY_TSV}"

conda run -n "${BIOM_ENV}" biom --help >/dev/null 2>&1 \
    || die "La commande biom est absente de l'environnement ${BIOM_ENV}."

log "Conversion de la table BIOM en TSV"
conda run -n "${BIOM_ENV}" biom convert \
    -i "${BIOM_TABLE}" \
    -o "${TSV_TABLE}" \
    --to-tsv

# La premiere ligne est un commentaire biom. La deuxieme est l'en-tete.
sed '1d; s/^#OTU ID/ASV_ID/' "${TSV_TABLE}" > "${EUK_ASV_TSV}"
[[ -s "${EUK_ASV_TSV}" ]] || die "Table ASV TSV absente ou vide : ${EUK_ASV_TSV}"

log "Fusion de la table ASV et de la taxonomie PR2"
python3 - "${EUK_ASV_TSV}" "${TAXONOMY_TSV}" "${EUK_ASV_TAXONOMY_TSV}" <<'PY'
import csv
import sys

asv_path, taxonomy_path, output_path = sys.argv[1:]

taxonomy = {}
confidence = {}
with open(taxonomy_path, encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        taxonomy[row["Feature ID"]] = row.get("Taxon", "Unassigned")
        confidence[row["Feature ID"]] = row.get("Confidence", "")

with open(asv_path, encoding="utf-8", newline="") as source, open(output_path, "w", encoding="utf-8", newline="") as destination:
    reader = csv.reader(source, delimiter="\t")
    writer = csv.writer(destination, delimiter="\t", lineterminator="\n")

    header = next(reader)
    writer.writerow(["ASV_ID", "PR2_taxonomy", "Confidence"] + header[1:])

    for row in reader:
        asv_id = row[0]
        writer.writerow([
            asv_id,
            taxonomy.get(asv_id, "Unassigned"),
            confidence.get(asv_id, ""),
        ] + row[1:])
PY

[[ -s "${EUK_ASV_TAXONOMY_TSV}" ]] \
    || die "Fusion ASV/taxonomie non produite : ${EUK_ASV_TAXONOMY_TSV}"

# ==================== FIN ====================

log "Classification eucaryote terminee avec succes."
log "Table candidats eucaryotes : ${EUK_TABLE}"
log "Sequences candidates eucaryotes : ${EUK_REP_SEQS}"
log "Taxonomie PR2 : ${EUK_TAXONOMY}"
log "Visualisation table : ${EUK_TABLE_QZV}"
log "Visualisation taxonomie : ${EUK_TAXONOMY_QZV}"
log "Barplot eucaryotes : ${EUK_BAR_PLOT_QZV}"
log "Table eucaryotes TSV : ${EUK_ASV_TSV}"
log "Taxonomie PR2 TSV : ${TAXONOMY_TSV}"
log "Table ASV + taxonomie PR2 TSV : ${EUK_ASV_TAXONOMY_TSV}"
