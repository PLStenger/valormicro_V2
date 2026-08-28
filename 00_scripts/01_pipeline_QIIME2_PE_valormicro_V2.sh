#!/usr/bin/env bash
# =============================================================================
# Pipeline paired-end : FastQC/MultiQC -> Trimmomatic -> FastQC/MultiQC
# -> QIIME 2 import -> DADA2 -> filtrage ASV -> arbre MAFFT/FastTree.
#
# Les echantillons analyses sont EXCLUSIVEMENT ceux avec R1, R2 et New_label
# renseignes dans le fichier 00_infos_data.xlsx. Les lignes vides/commentaires
# du tableur sont ignorees. New_label devient l'identifiant QIIME 2.
#
# Lancement : nohup bash 01_pipeline_QIIME2_PE_valormicro_V2_corrige.sh \
#              > pipeline.log 2>&1 &
# =============================================================================

set -Eeuo pipefail
shopt -s nullglob
IFS=$'\n\t'

# Les scripts d'activation OpenJDK/Conda peuvent lire ces variables alors
# qu'elles ne sont pas predefinies. Les initialiser permet de conserver set -u.
export JAVA_HOME="${JAVA_HOME:-}"
export JAVA_LD_LIBRARY_PATH="${JAVA_LD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# Affiche la commande/l'identifiant de ligne en cas d'echec non gere.
trap 'rc=$?; echo "[ERREUR] Code ${rc}, ligne ${LINENO}: ${BASH_COMMAND}" >&2; exit "${rc}"' ERR

# ------------------------------- CONFIGURATION ------------------------------
PROJECT_DIR="/nvme/bio/data_fungi/valormicro_V2"
RAW_DIR="${PROJECT_DIR}/01_raw_data"
INFO_XLSX="${RAW_DIR}/00_infos_data.xlsx"
RESULTS_DIR="${PROJECT_DIR}/02_amplicon_pipeline"

THREADS=8
TRIMMOMATIC_HEAP="60G"
QIIME_THREADS=8
TMPDIR_BASE="${PROJECT_DIR}/tmp"

# Adapter si les noms de vos environnements Conda different.
EXCEL_ENV="excel_tools"
FASTQC_ENV="fastqc"
MULTIQC_ENV="multiqc"
TRIMMOMATIC_ENV="trimmomatic"
QIIME2_ENV="qiime2-2021.4"

# Fichier d'adaptateurs requis si TRIMMOMATIC_ADAPTERS=true.
ADAPTER_FILE="${PROJECT_DIR}/99_softwares/adapters_sequences.fasta"
TRIMMOMATIC_ADAPTERS=true

# Parametres Trimmomatic, repris de votre ancien pipeline.
LEADING=30
TRAILING=30
SLIDINGWINDOW="26:30"
MINLEN=150

# Parametres DADA2. A determiner apres inspection de demux.qzv. Une valeur 0
# pour --p-trunc-len-* conserve la longueur complete.
DADA2_TRIM_LEFT_F=0
DADA2_TRIM_LEFT_R=0
DADA2_TRUNC_LEN_F=0
DADA2_TRUNC_LEN_R=0
DADA2_MAX_EE_F=2
DADA2_MAX_EE_R=2
DADA2_CHIM_METHOD="consensus"
MIN_SAMPLES_FEATURE=2

# Active/desactive les grandes etapes.
RUN_FASTQC_RAW=true
RUN_TRIMMOMATIC=true
RUN_FASTQC_CLEAN=true
RUN_QIIME_IMPORT=true
RUN_DADA2=true
RUN_FILTER=true
RUN_TREE=true

# --------------------------------- CHEMINS -----------------------------------
QC_RAW_DIR="${RESULTS_DIR}/01_qc_raw"
CLEAN_DIR="${RESULTS_DIR}/02_cleaned_data"
QC_CLEAN_DIR="${RESULTS_DIR}/03_qc_cleaned"
DATABASE_DIR="${RESULTS_DIR}/04_database_files"
QIIME_DIR="${RESULTS_DIR}/05_qiime2"
QIIME_CORE="${QIIME_DIR}/core"
QIIME_VISUAL="${QIIME_DIR}/visual"
QIIME_TREE="${QIIME_DIR}/tree"
QIIME_EXPORT="${QIIME_DIR}/export"
LOG_DIR="${RESULTS_DIR}/logs"

SAMPLE_SHEET="${DATABASE_DIR}/samples.tsv"
RAW_FASTQ_LIST="${DATABASE_DIR}/raw_fastq_list.txt"
RAW_PAIRS_TSV="${DATABASE_DIR}/raw_pairs.tsv"
MANIFEST="${DATABASE_DIR}/manifest_pe.tsv"
METADATA="${DATABASE_DIR}/sample-metadata.tsv"

log() {
    printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOG_DIR}/pipeline.log"
}

die() {
    log "ERREUR : $*"
    exit 1
}

activate_env() {
    local env_name="$1"
    # Protection a chaque conda activate, y compris quand les activateurs
    # externes OpenJDK sont evalues sous set -u.
    export JAVA_HOME="${JAVA_HOME:-}"
    export JAVA_LD_LIBRARY_PATH="${JAVA_LD_LIBRARY_PATH:-}"
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
    conda deactivate >/dev/null 2>&1 || true
    conda activate "${env_name}"
}

check_command() {
    command -v "$1" >/dev/null 2>&1 || die "Commande introuvable : $1"
}

run_fastqc_multiqc_from_list() {
    local input_list="$1"
    local output_dir="$2"
    local label="$3"
    local files=()
    mapfile -t files < "${input_list}"
    ((${#files[@]} > 0)) || die "La liste FastQC est vide : ${input_list}"

    mkdir -p "${output_dir}/fastqc" "${output_dir}/multiqc"
    log "FastQC (${label}) sur ${#files[@]} fichiers declares dans le tableur"

    activate_env "${FASTQC_ENV}"
    check_command fastqc
    fastqc --threads "${THREADS}" --outdir "${output_dir}/fastqc" "${files[@]}"

    activate_env "${MULTIQC_ENV}"
    check_command multiqc
    multiqc --force --outdir "${output_dir}/multiqc" "${output_dir}/fastqc"
}

run_fastqc_multiqc_from_directory() {
    local input_dir="$1"
    local output_dir="$2"
    local label="$3"
    local files=("${input_dir}"/*.fastq.gz)
    ((${#files[@]} > 0)) || die "Aucun fichier *.fastq.gz dans ${input_dir}"

    mkdir -p "${output_dir}/fastqc" "${output_dir}/multiqc"
    log "FastQC (${label}) sur ${#files[@]} fichiers"

    activate_env "${FASTQC_ENV}"
    check_command fastqc
    fastqc --threads "${THREADS}" --outdir "${output_dir}/fastqc" "${files[@]}"

    activate_env "${MULTIQC_ENV}"
    check_command multiqc
    multiqc --force --outdir "${output_dir}/multiqc" "${output_dir}/fastqc"
}

CONDA_BASE="$(conda info --base 2>/dev/null || true)"
if [[ -z "${CONDA_BASE}" || ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    echo "ERREUR : Conda est introuvable. Chargez miniconda/anaconda avant de lancer ce script." >&2
    exit 1
fi
source "${CONDA_BASE}/etc/profile.d/conda.sh"

mkdir -p "${QC_RAW_DIR}" "${CLEAN_DIR}" "${QC_CLEAN_DIR}" "${DATABASE_DIR}" \
    "${QIIME_CORE}" "${QIIME_VISUAL}" "${QIIME_TREE}" "${QIIME_EXPORT}" \
    "${LOG_DIR}" "${TMPDIR_BASE}"
export TMPDIR="${TMPDIR_BASE}"

[[ -f "${INFO_XLSX}" ]] || die "Tableur absent : ${INFO_XLSX}"
log "Demarrage du pipeline dans ${RESULTS_DIR}"
log "Fichier d'informations : ${INFO_XLSX}"

# 1. Lire le tableur une seule fois et creer toutes les tables necessaires.
log "Creation de samples.tsv, raw_fastq_list.txt, raw_pairs.tsv, manifest_pe.tsv et sample-metadata.tsv"
activate_env "${EXCEL_ENV}"
check_command python

python - "${INFO_XLSX}" "${SAMPLE_SHEET}" "${RAW_FASTQ_LIST}" "${RAW_PAIRS_TSV}" "${MANIFEST}" "${METADATA}" "${RAW_DIR}" "${CLEAN_DIR}" <<'PY'
import re
import sys
from pathlib import Path

import openpyxl  # Verification explicite de la dependance du lecteur xlsx.
import pandas as pd

xlsx, sample_sheet, raw_fastq_list, raw_pairs_tsv, manifest, metadata, raw_dir, clean_dir = sys.argv[1:]
df = pd.read_excel(xlsx, dtype=str)
df.columns = [str(col).strip() for col in df.columns]
required = ["R1", "R2", "New_label"]
missing = [column for column in required if column not in df.columns]
if missing:
    raise SystemExit(
        "Colonnes absentes du tableur : " + ", ".join(missing) +
        ". Colonnes disponibles : " + ", ".join(df.columns)
    )

# Une ligne est un echantillon uniquement si les trois champs obligatoires sont
# complets. Les lignes totalement vides ou de commentaire sont ignorees.
for column in required:
    df[column] = df[column].fillna("").astype(str).str.strip()

has_any_required_value = (df["R1"] != "") | (df["R2"] != "") | (df["New_label"] != "")
partially_filled = df.loc[
    has_any_required_value &
    ((df["R1"] == "") | (df["R2"] == "") | (df["New_label"] == "")),
    ["R1", "R2", "New_label"]
]
if not partially_filled.empty:
    raise SystemExit(
        "Lignes partiellement renseignees : R1, R2 et New_label sont obligatoires.\n" +
        partially_filled.to_string()
    )

df = df.loc[(df["R1"] != "") & (df["R2"] != "") & (df["New_label"] != "")].copy()
if df.empty:
    raise SystemExit("Aucun echantillon valide (R1/R2/New_label) dans le tableur.")

# Les IDs QIIME 2 sont issus de New_label, avec remplacement des espaces par _.
df["sample-id"] = df["New_label"].str.replace(r"\s+", "_", regex=True)
valid_id = r"^[A-Za-z0-9_.-]+$"
bad_ids = df.loc[~df["sample-id"].str.match(valid_id), "sample-id"].tolist()
if bad_ids:
    raise SystemExit("Identifiants QIIME 2 non valides : " + ", ".join(bad_ids))
if df["sample-id"].duplicated().any():
    duplicated = df.loc[df["sample-id"].duplicated(keep=False), "sample-id"].tolist()
    raise SystemExit("New_label dupliques : " + ", ".join(duplicated))

# Table d'audit complete.
df.to_csv(sample_sheet, sep="\t", index=False, na_rep="")

raw_pairs = pd.DataFrame({
    "sample-id": df["sample-id"],
    "raw-forward-absolute-filepath": [str(Path(raw_dir) / name) for name in df["R1"]],
    "raw-reverse-absolute-filepath": [str(Path(raw_dir) / name) for name in df["R2"]],
})
raw_pairs.to_csv(raw_pairs_tsv, sep="\t", index=False)

# Cette liste ne contient QUE les 2 FASTQ declares par echantillon : FastQC raw
# ne scannera plus les autres runs ou projets presents dans RAW_DIR.
with open(raw_fastq_list, "w", encoding="utf-8") as handle:
    for path in raw_pairs["raw-forward-absolute-filepath"]:
        handle.write(path + "\n")
    for path in raw_pairs["raw-reverse-absolute-filepath"]:
        handle.write(path + "\n")

# Manifest QIIME 2 : reads paired apres Trimmomatic.
manifest_df = pd.DataFrame({
    "sample-id": df["sample-id"],
    "forward-absolute-filepath": [str(Path(clean_dir) / f"{sid}_R1_paired.fastq.gz") for sid in df["sample-id"]],
    "reverse-absolute-filepath": [str(Path(clean_dir) / f"{sid}_R2_paired.fastq.gz") for sid in df["sample-id"]],
})
manifest_df.to_csv(manifest, sep="\t", index=False)

# Metadata compatible QIIME 2 : conserver les informations descriptives, mais
# exclure les noms de FASTQ et New_label deja represente par #SampleID.
metadata_df = df.drop(columns=["R1", "R2", "New_label", "sample-id"], errors="ignore").copy()
metadata_df.insert(0, "#SampleID", df["sample-id"])
metadata_df.columns = [
    "#SampleID" if column == "#SampleID" else re.sub(r"#", "", str(column)).strip().replace(" ", "_")
    for column in metadata_df.columns
]
metadata_df.to_csv(metadata, sep="\t", index=False, na_rep="")

print(f"{len(df)} echantillons ecrits")
PY

[[ -s "${SAMPLE_SHEET}" ]] || die "samples.tsv absent ou vide apres lecture Excel : ${SAMPLE_SHEET}"
[[ -s "${RAW_FASTQ_LIST}" ]] || die "raw_fastq_list.txt absent ou vide : ${RAW_FASTQ_LIST}"
[[ -s "${RAW_PAIRS_TSV}" ]] || die "raw_pairs.tsv absent ou vide : ${RAW_PAIRS_TSV}"
[[ -s "${MANIFEST}" ]] || die "manifest_pe.tsv absent ou vide : ${MANIFEST}"
[[ -s "${METADATA}" ]] || die "sample-metadata.tsv absent ou vide : ${METADATA}"

N_SAMPLES=$(( $(wc -l < "${MANIFEST}") - 1 ))
N_RAW_FASTQ=$(wc -l < "${RAW_FASTQ_LIST}")
(( N_SAMPLES > 0 )) || die "Le manifest ne contient aucun echantillon."
(( N_RAW_FASTQ == N_SAMPLES * 2 )) || die "La liste FASTQ ne contient pas exactement deux fichiers par echantillon."
log "${N_SAMPLES} echantillons trouves dans le tableur ; ${N_RAW_FASTQ} FASTQ bruts seront analyses"

# 2. Controle explicite de la presence et de l'integrite de chaque FASTQ brut.
log "Validation des FASTQ bruts declares"
while IFS=$'\t' read -r sample_id raw_r1 raw_r2; do
    [[ "${sample_id}" == "sample-id" ]] && continue
    [[ -s "${raw_r1}" ]] || die "R1 introuvable ou vide pour ${sample_id} : ${raw_r1}"
    [[ -s "${raw_r2}" ]] || die "R2 introuvable ou vide pour ${sample_id} : ${raw_r2}"
    gzip -t "${raw_r1}"
    gzip -t "${raw_r2}"
done < "${RAW_PAIRS_TSV}"

# 3. FastQC/MultiQC uniquement sur les reads R1/R2 declares dans le tableur.
if [[ "${RUN_FASTQC_RAW}" == true ]]; then
    run_fastqc_multiqc_from_list "${RAW_FASTQ_LIST}" "${QC_RAW_DIR}" "reads bruts declares"
fi

# 4. Trimmomatic PE. Cette boucle lit raw_pairs.tsv et ne depend pas de pandas :
# l'environnement trimmomatic n'a donc pas besoin de pandas.
if [[ "${RUN_TRIMMOMATIC}" == true ]]; then
    if [[ "${TRIMMOMATIC_ADAPTERS}" == true ]]; then
        [[ -f "${ADAPTER_FILE}" ]] || die "Fichier d'adaptateurs absent : ${ADAPTER_FILE}"
        TRIM_ADAPTER_ARG=("ILLUMINACLIP:${ADAPTER_FILE}:2:30:10")
    else
        TRIM_ADAPTER_ARG=()
    fi

    activate_env "${TRIMMOMATIC_ENV}"
    check_command trimmomatic
    log "Trimmomatic PE sur ${N_SAMPLES} echantillons"

    while IFS=$'\t' read -r sample_id raw_r1 raw_r2; do
        [[ "${sample_id}" == "sample-id" ]] && continue

        out_r1_paired="${CLEAN_DIR}/${sample_id}_R1_paired.fastq.gz"
        out_r1_unpaired="${CLEAN_DIR}/${sample_id}_R1_unpaired.fastq.gz"
        out_r2_paired="${CLEAN_DIR}/${sample_id}_R2_paired.fastq.gz"
        out_r2_unpaired="${CLEAN_DIR}/${sample_id}_R2_unpaired.fastq.gz"

        if [[ -s "${out_r1_paired}" && -s "${out_r2_paired}" ]]; then
            log "Trimmomatic deja termine : ${sample_id}"
            continue
        fi

        log "Trimmomatic : ${sample_id}"
        trimmomatic PE -Xmx"${TRIMMOMATIC_HEAP}" -threads "${THREADS}" -phred33 \
            "${raw_r1}" "${raw_r2}" \
            "${out_r1_paired}" "${out_r1_unpaired}" \
            "${out_r2_paired}" "${out_r2_unpaired}" \
            "${TRIM_ADAPTER_ARG[@]}" \
            "LEADING:${LEADING}" \
            "TRAILING:${TRAILING}" \
            "SLIDINGWINDOW:${SLIDINGWINDOW}" \
            "MINLEN:${MINLEN}"
    done < "${RAW_PAIRS_TSV}"
fi

# Verifie que les reads paired necessaires au manifest existent apres trimming.
while IFS=$'\t' read -r sample_id paired_r1 paired_r2; do
    [[ "${sample_id}" == "sample-id" ]] && continue
    [[ -s "${paired_r1}" ]] || die "R1 paired manquant pour ${sample_id} : ${paired_r1}"
    [[ -s "${paired_r2}" ]] || die "R2 paired manquant pour ${sample_id} : ${paired_r2}"
done < "${MANIFEST}"

# 5. QC post-trimming. Ici tous les *.fastq.gz de CLEAN_DIR proviennent de ce
# pipeline : paired ET unpaired sont inclus afin d'evaluer le trimming.
if [[ "${RUN_FASTQC_CLEAN}" == true ]]; then
    run_fastqc_multiqc_from_directory "${CLEAN_DIR}" "${QC_CLEAN_DIR}" "reads nettoyes"
fi

# 6. Import QIIME 2.
activate_env "${QIIME2_ENV}"
check_command qiime
export TMPDIR="${TMPDIR_BASE}"

if [[ "${RUN_QIIME_IMPORT}" == true ]]; then
    log "Import QIIME 2 du manifest paired-end"
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "${MANIFEST}" \
        --input-format PairedEndFastqManifestPhred33V2 \
        --output-path "${QIIME_CORE}/demux.qza"

    qiime demux summarize \
        --i-data "${QIIME_CORE}/demux.qza" \
        --o-visualization "${QIIME_VISUAL}/demux.qzv"
fi

[[ -f "${QIIME_CORE}/demux.qza" ]] || die "demux.qza absent : activez RUN_QIIME_IMPORT ou fournissez ce fichier."

# 7. Denoising DADA2.
if [[ "${RUN_DADA2}" == true ]]; then
    log "Denoising DADA2 paired-end"
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "${QIIME_CORE}/demux.qza" \
        --p-trim-left-f "${DADA2_TRIM_LEFT_F}" \
        --p-trim-left-r "${DADA2_TRIM_LEFT_R}" \
        --p-trunc-len-f "${DADA2_TRUNC_LEN_F}" \
        --p-trunc-len-r "${DADA2_TRUNC_LEN_R}" \
        --p-max-ee-f "${DADA2_MAX_EE_F}" \
        --p-max-ee-r "${DADA2_MAX_EE_R}" \
        --p-chimera-method "${DADA2_CHIM_METHOD}" \
        --p-n-threads "${QIIME_THREADS}" \
        --o-table "${QIIME_CORE}/table.qza" \
        --o-representative-sequences "${QIIME_CORE}/rep-seqs.qza" \
        --o-denoising-stats "${QIIME_CORE}/denoising-stats.qza"

    qiime metadata tabulate \
        --m-input-file "${QIIME_CORE}/denoising-stats.qza" \
        --o-visualization "${QIIME_VISUAL}/denoising-stats.qzv"

    qiime feature-table summarize \
        --i-table "${QIIME_CORE}/table.qza" \
        --m-sample-metadata-file "${METADATA}" \
        --o-visualization "${QIIME_VISUAL}/table.qzv"

    qiime feature-table tabulate-seqs \
        --i-data "${QIIME_CORE}/rep-seqs.qza" \
        --o-visualization "${QIIME_VISUAL}/rep-seqs.qzv"
fi

[[ -f "${QIIME_CORE}/table.qza" ]] || die "table.qza absent : activez RUN_DADA2 ou fournissez ce fichier."
[[ -f "${QIIME_CORE}/rep-seqs.qza" ]] || die "rep-seqs.qza absent : activez RUN_DADA2 ou fournissez ce fichier."

# 8. Filtrage de contingence. Aucun filtrage decontam automatique n'est realise
# car aucun controle negatif n'est identifie dans le tableur.
if [[ "${RUN_FILTER}" == true ]]; then
    log "Filtrage des ASV presentes dans au moins ${MIN_SAMPLES_FEATURE} echantillons"
    qiime feature-table filter-features \
        --i-table "${QIIME_CORE}/table.qza" \
        --p-min-samples "${MIN_SAMPLES_FEATURE}" \
        --o-filtered-table "${QIIME_CORE}/table_min${MIN_SAMPLES_FEATURE}samples.qza"

    qiime feature-table filter-seqs \
        --i-data "${QIIME_CORE}/rep-seqs.qza" \
        --i-table "${QIIME_CORE}/table_min${MIN_SAMPLES_FEATURE}samples.qza" \
        --o-filtered-data "${QIIME_CORE}/rep-seqs_min${MIN_SAMPLES_FEATURE}samples.qza"

    qiime feature-table summarize \
        --i-table "${QIIME_CORE}/table_min${MIN_SAMPLES_FEATURE}samples.qza" \
        --m-sample-metadata-file "${METADATA}" \
        --o-visualization "${QIIME_VISUAL}/table_min${MIN_SAMPLES_FEATURE}samples.qzv"
fi

FILTERED_REPSEQS="${QIIME_CORE}/rep-seqs_min${MIN_SAMPLES_FEATURE}samples.qza"
[[ -f "${FILTERED_REPSEQS}" ]] || FILTERED_REPSEQS="${QIIME_CORE}/rep-seqs.qza"

# 9. Arbre ASV : MAFFT -> mask -> FastTree -> midpoint root.
if [[ "${RUN_TREE}" == true ]]; then
    log "Construction de l'arbre phylogenetique des ASV"
    qiime alignment mafft \
        --i-sequences "${FILTERED_REPSEQS}" \
        --p-n-threads "${QIIME_THREADS}" \
        --o-alignment "${QIIME_TREE}/aligned-rep-seqs.qza"

    qiime alignment mask \
        --i-alignment "${QIIME_TREE}/aligned-rep-seqs.qza" \
        --o-masked-alignment "${QIIME_TREE}/masked-aligned-rep-seqs.qza"

    qiime phylogeny fasttree \
        --i-alignment "${QIIME_TREE}/masked-aligned-rep-seqs.qza" \
        --o-tree "${QIIME_TREE}/unrooted-tree.qza"

    qiime phylogeny midpoint-root \
        --i-tree "${QIIME_TREE}/unrooted-tree.qza" \
        --o-rooted-tree "${QIIME_TREE}/rooted-tree.qza"
fi

# 10. Exports QIIME 2.
log "Export des resultats QIIME 2"
mkdir -p "${QIIME_EXPORT}/core" "${QIIME_EXPORT}/tree" "${QIIME_EXPORT}/visual"

qiime tools export --input-path "${QIIME_CORE}/table.qza" --output-path "${QIIME_EXPORT}/core/table"
qiime tools export --input-path "${QIIME_CORE}/rep-seqs.qza" --output-path "${QIIME_EXPORT}/core/rep-seqs"
qiime tools export --input-path "${QIIME_CORE}/denoising-stats.qza" --output-path "${QIIME_EXPORT}/core/denoising-stats"

if [[ -f "${QIIME_CORE}/table_min${MIN_SAMPLES_FEATURE}samples.qza" ]]; then
    qiime tools export \
        --input-path "${QIIME_CORE}/table_min${MIN_SAMPLES_FEATURE}samples.qza" \
        --output-path "${QIIME_EXPORT}/core/table_min${MIN_SAMPLES_FEATURE}samples"
    qiime tools export \
        --input-path "${QIIME_CORE}/rep-seqs_min${MIN_SAMPLES_FEATURE}samples.qza" \
        --output-path "${QIIME_EXPORT}/core/rep-seqs_min${MIN_SAMPLES_FEATURE}samples"
fi

if [[ -f "${QIIME_TREE}/rooted-tree.qza" ]]; then
    qiime tools export --input-path "${QIIME_TREE}/aligned-rep-seqs.qza" --output-path "${QIIME_EXPORT}/tree/aligned-rep-seqs"
    qiime tools export --input-path "${QIIME_TREE}/masked-aligned-rep-seqs.qza" --output-path "${QIIME_EXPORT}/tree/masked-aligned-rep-seqs"
    qiime tools export --input-path "${QIIME_TREE}/unrooted-tree.qza" --output-path "${QIIME_EXPORT}/tree/unrooted-tree"
    qiime tools export --input-path "${QIIME_TREE}/rooted-tree.qza" --output-path "${QIIME_EXPORT}/tree/rooted-tree"
fi

for qzv in "${QIIME_VISUAL}"/*.qzv; do
    [[ -e "${qzv}" ]] || continue
    qzv_name="$(basename "${qzv}" .qzv)"
    qiime tools export --input-path "${qzv}" --output-path "${QIIME_EXPORT}/visual/${qzv_name}"
done

log "Pipeline termine avec succes."
log "Manifest : ${MANIFEST}"
log "Metadata : ${METADATA}"
log "Resultats QIIME 2 : ${QIIME_DIR}"
