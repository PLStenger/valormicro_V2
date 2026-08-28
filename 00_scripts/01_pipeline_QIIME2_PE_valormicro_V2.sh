#!/usr/bin/env bash
# ============================================================================
# Pipeline PE : FastQC/MultiQC -> Trimmomatic -> FastQC/MultiQC -> QIIME 2
# -> DADA2 -> arbre phylogenetique.
# Projet : valormicro_V2
#
# Les manifestes et le fichier metadata sont construits automatiquement depuis
# 00_infos_data.xlsx. Les identifiants QIIME 2 sont les valeurs de New_label.
#
# Prerequis : bash, awk, gzip, Python >= 3.8 avec pandas + openpyxl,
# FastQC, MultiQC, Trimmomatic et QIIME 2 dans les environnements Conda definis
# ci-dessous. Adapter les variables de CONFIGURATION avant execution.
#
# Lancer : bash 01_pipeline_QIIME2_PE_valormicro_V2.sh
# Reprendre une execution : bash 01_pipeline_QIIME2_PE_valormicro_V2.sh
# ============================================================================

IFS=$'\n\t'

# -------------------------------- CONFIGURATION -----------------------------
PROJECT_DIR="/nvme/bio/data_fungi/valormicro_V2"
RAW_DIR="${PROJECT_DIR}/01_raw_data"
INFO_XLSX="${RAW_DIR}/00_infos_data.xlsx"
RESULTS_DIR="${PROJECT_DIR}/02_amplicon_pipeline"

# Ressources.
THREADS=8
TRIMMOMATIC_HEAP="60G"
QIIME_THREADS=8
TMPDIR_BASE="${PROJECT_DIR}/tmp"

# Environnements Conda : modifier ces noms si necessaire.
FASTQC_ENV="fastqc"
MULTIQC_ENV="multiqc"
TRIMMOMATIC_ENV="trimmomatic"
QIIME2_ENV="qiime2-2021.4"

# Fichier FASTA d'adaptateurs Trimmomatic. Obligatoire seulement si
# TRIMMOMATIC_ADAPTERS=true.
ADAPTER_FILE="${PROJECT_DIR}/99_softwares/adapters_sequences.fasta"
TRIMMOMATIC_ADAPTERS=true

# Parametres de trimming (reprennent les anciens scripts).
LEADING=30
TRAILING=30
SLIDINGWINDOW="26:30"
MINLEN=150

# Parametres DADA2. Par defaut, aucun rognage supplementaire : choisir les
# valeurs a partir de 03_qc_cleaned/multiqc_report.html ou demux.qzv.
# Mettre 0 pour --p-trunc-len-* afin de conserver la longueur complete.
DADA2_TRIM_LEFT_F=0
DADA2_TRIM_LEFT_R=0
DADA2_TRUNC_LEN_F=0
DADA2_TRUNC_LEN_R=0
DADA2_MAX_EE_F=2
DADA2_MAX_EE_R=2
DADA2_CHIM_METHOD="consensus"

# Active desactive les etapes. L'arbre est realise sur les ASV filtrees pour
# occurrence dans >= MIN_SAMPLES_FEATURE echantillons.
RUN_FASTQC_RAW=true
RUN_TRIMMOMATIC=true
RUN_FASTQC_CLEAN=true
RUN_QIIME_IMPORT=true
RUN_DADA2=true
RUN_FILTER=true
RUN_TREE=true
MIN_SAMPLES_FEATURE=2

# ----------------------------------------------------------------------------
QC_RAW_DIR="${RESULTS_DIR}/01_qc_raw"
CLEAN_DIR="${RESULTS_DIR}/02_cleaned_data"
QC_CLEAN_DIR="${RESULTS_DIR}/03_qc_cleaned"
DATABASE_DIR="${RESULTS_DIR}/04_database_files"
QIIME_DIR="${RESULTS_DIR}/05_qiime2"
QIIME_CORE="${QIIME_DIR}/core"
QIIME_VISUAL="${QIIME_DIR}/visual"
QIIME_TREE="${QIIME_DIR}/tree"
QIIME_EXPORT="${QIIME_DIR}/export"
MANIFEST="${DATABASE_DIR}/manifest_pe.tsv"
METADATA="${DATABASE_DIR}/sample-metadata.tsv"
SAMPLE_SHEET="${DATABASE_DIR}/samples.tsv"
LOG_DIR="${RESULTS_DIR}/logs"

CONDA_BASE="$(conda info --base 2>/dev/null || true)"
if [[ -z "${CONDA_BASE}" || ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    echo "ERREUR : Conda est introuvable. Chargez/miniconda avant de lancer ce script." >&2
    exit 1
fi
source "${CONDA_BASE}/etc/profile.d/conda.sh"

log() {
    printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOG_DIR}/pipeline.log"
}

die() {
    log "ERREUR : $*"
    exit 1
}

activate_env() {
    conda deactivate >/dev/null 2>&1 || true
    conda activate "$1"
}

check_command() {
    command -v "$1" >/dev/null 2>&1 || die "Commande introuvable : $1"
}

run_fastqc_multiqc() {
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

mkdir -p "${QC_RAW_DIR}" "${CLEAN_DIR}" "${QC_CLEAN_DIR}" "${DATABASE_DIR}" \
         "${QIIME_CORE}" "${QIIME_VISUAL}" "${QIIME_TREE}" "${QIIME_EXPORT}" \
         "${LOG_DIR}" "${TMPDIR_BASE}"
export TMPDIR="${TMPDIR_BASE}"

[[ -f "${INFO_XLSX}" ]] || die "Tableur absent : ${INFO_XLSX}"
log "Demarrage du pipeline dans ${RESULTS_DIR}"
log "Fichier d'informations : ${INFO_XLSX}"

# 1. Conversion de l'Excel en TSV, manifeste PE et metadata QIIME 2.
log "Creation de samples.tsv, manifest_pe.tsv et sample-metadata.tsv"
python3 - "${INFO_XLSX}" "${SAMPLE_SHEET}" "${MANIFEST}" "${METADATA}" "${CLEAN_DIR}" <<'PY'
import os
import re
import sys
import unicodedata
from pathlib import Path

import pandas as pd

xlsx, sample_sheet, manifest, metadata, clean_dir = sys.argv[1:]
df = pd.read_excel(xlsx, dtype=str)
df.columns = [str(c).strip() for c in df.columns]

required = ["R1", "R2", "New_label"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit(f"Colonnes absentes du tableur : {', '.join(missing)}. Colonnes disponibles : {', '.join(df.columns)}")

# Enleve les lignes entierement vides puis verifie les champs indispensables.
df = df.dropna(how="all").copy()
for col in required:
    df[col] = df[col].fillna("").astype(str).str.strip()
invalid = df.index[(df["R1"] == "") | (df["R2"] == "") | (df["New_label"] == "")].tolist()
if invalid:
    raise SystemExit(f"Lignes Excel incompletes (index pandas) : {invalid}")

# Les ID QIIME 2 doivent etre uniques et ne pas contenir de blancs.
df["sample-id"] = df["New_label"].str.replace(r"\s+", "_", regex=True)
if df["sample-id"].duplicated().any():
    dup = df.loc[df["sample-id"].duplicated(keep=False), "sample-id"].tolist()
    raise SystemExit(f"New_label dupliques : {dup}")
if df["sample-id"].str.contains(r"[^A-Za-z0-9_.-]", regex=True).any():
    bad = df.loc[df["sample-id"].str.contains(r"[^A-Za-z0-9_.-]", regex=True), "sample-id"].tolist()
    raise SystemExit(f"ID QIIME 2 non valides apres normalisation : {bad}")

# Table de correspondance complete (audit/reproductibilite).
df.to_csv(sample_sheet, sep="\t", index=False, na_rep="")

# Le manifest reference les reads pairs produits par Trimmomatic.
manifest_df = pd.DataFrame({
    "sample-id": df["sample-id"],
    "forward-absolute-filepath": [str(Path(clean_dir) / f"{sid}_R1_paired.fastq.gz") for sid in df["sample-id"]],
    "reverse-absolute-filepath": [str(Path(clean_dir) / f"{sid}_R2_paired.fastq.gz") for sid in df["sample-id"]],
})
manifest_df.to_csv(manifest, sep="\t", index=False)

# Metadata QIIME 2 : toutes les colonnes du tableur sauf les noms de fichiers.
meta = df.drop(columns=["R1", "R2", "New_label"], errors="ignore").copy()
meta.insert(0, "#SampleID", df["sample-id"])
meta = meta.drop(columns=["sample-id"], errors="ignore")
# Les caracteres '#' sont reserves dans un fichier metadata QIIME 2.
meta.columns = ["#SampleID"] + [re.sub(r"#", "", str(c)).strip().replace(" ", "_") for c in meta.columns[1:]]
meta.to_csv(metadata, sep="\t", index=False, na_rep="")

print(f"{len(df)} echantillons ecrits")
PY

N_SAMPLES=$(( $(wc -l < "${MANIFEST}") - 1 ))
((N_SAMPLES > 0)) || die "Le manifest ne contient aucun echantillon."
log "${N_SAMPLES} echantillons trouves dans le tableur"

# Verifie la presence et l'integrite gzip des FASTQ bruts declares dans l'Excel.
log "Validation des FASTQ bruts declares"
while IFS=$'\t' read -r sample_id r1 r2; do
    [[ "${sample_id}" == "sample-id" ]] && continue
    # Les chemins sources sont recuperes de samples.tsv, pas du manifest.
    :
done < "${MANIFEST}"

while IFS=$'\t' read -r sample_id r1 r2; do
    [[ "${sample_id}" == "sample-id" ]] && continue
    # Aucune action dans cette boucle : elle conserve un controle explicite du manifest.
    [[ -n "${r1}" && -n "${r2}" ]] || die "Manifest incomplet pour ${sample_id}"
done < "${MANIFEST}"

# Utilise samples.tsv pour relier les fichiers bruts aux identifiants New_label.
while IFS=$'\t' read -r label site gps r1 r2 cluster_path new_label rest; do
    [[ "${label}" == "Label" ]] && continue
    [[ -n "${r1}" && -n "${r2}" && -n "${new_label}" ]] || continue
    raw_r1="${RAW_DIR}/${r1}"
    raw_r2="${RAW_DIR}/${r2}"
    [[ -s "${raw_r1}" ]] || die "R1 introuvable ou vide : ${raw_r1}"
    [[ -s "${raw_r2}" ]] || die "R2 introuvable ou vide : ${raw_r2}"
    gzip -t "${raw_r1}"
    gzip -t "${raw_r2}"
done < <(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++){h[$i]=i}; next} {print $(h["Label"]),$(h["Nom du lieu"]),$(h["Coordonnées GPS Google Maps"]),$(h["R1"]),$(h["R2"]),$(h["Cluster_path"]),$(h["New_label"]),$(h["R color (see https://usmap.dev/docs/Rcolor.pdf)"])}' "${SAMPLE_SHEET}")

if [[ "${RUN_FASTQC_RAW}" == true ]]; then
    run_fastqc_multiqc "${RAW_DIR}" "${QC_RAW_DIR}" "reads bruts"
fi

# 2. Trimming PE et renommage deterministe en identifiants New_label.
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

    python3 - "${SAMPLE_SHEET}" <<'PY' | while IFS=$'\t' read -r sample_id r1 r2; do
import sys
import pandas as pd

df = pd.read_csv(sys.argv[1], sep='\t', dtype=str).fillna('')
for _, row in df.iterrows():
    print(f"{row['sample-id']}\t{row['R1']}\t{row['R2']}")
PY
        raw_r1="${RAW_DIR}/${r1}"
        raw_r2="${RAW_DIR}/${r2}"
        out_r1p="${CLEAN_DIR}/${sample_id}_R1_paired.fastq.gz"
        out_r1u="${CLEAN_DIR}/${sample_id}_R1_unpaired.fastq.gz"
        out_r2p="${CLEAN_DIR}/${sample_id}_R2_paired.fastq.gz"
        out_r2u="${CLEAN_DIR}/${sample_id}_R2_unpaired.fastq.gz"

        if [[ -s "${out_r1p}" && -s "${out_r2p}" ]]; then
            log "Trimmomatic deja termine : ${sample_id}"
            continue
        fi

        log "Trimmomatic : ${sample_id}"
        trimmomatic PE -Xmx"${TRIMMOMATIC_HEAP}" -threads "${THREADS}" -phred33 \
            "${raw_r1}" "${raw_r2}" \
            "${out_r1p}" "${out_r1u}" "${out_r2p}" "${out_r2u}" \
            "${TRIM_ADAPTER_ARG[@]}" \
            "LEADING:${LEADING}" "TRAILING:${TRAILING}" \
            "SLIDINGWINDOW:${SLIDINGWINDOW}" "MINLEN:${MINLEN}"
    done
fi

# Les fichiers paires sont indispensables pour DADA2 PE.
while IFS=$'\t' read -r sample_id r1 r2; do
    [[ "${sample_id}" == "sample-id" ]] && continue
    [[ -s "${r1}" ]] || die "Read R1 paire manquant : ${r1}"
    [[ -s "${r2}" ]] || die "Read R2 paire manquant : ${r2}"
done < "${MANIFEST}"

if [[ "${RUN_FASTQC_CLEAN}" == true ]]; then
    run_fastqc_multiqc "${CLEAN_DIR}" "${QC_CLEAN_DIR}" "reads nettoyes"
fi

activate_env "${QIIME2_ENV}"
check_command qiime
export TMPDIR="${TMPDIR_BASE}"

# 3. Import QIIME 2 et resume de la qualite.
if [[ "${RUN_QIIME_IMPORT}" == true ]]; then
    log "Import QIIME 2 du manifest PE"
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

# 4. Denoising DADA2.
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

# 5. Filtrage de contingence. Aucune soustraction par controles negatifs n'est
# appliquee automatiquement car aucun controle negatif n'est identifie dans le
# tableur fourni.
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

FILTERED_TABLE="${QIIME_CORE}/table_min${MIN_SAMPLES_FEATURE}samples.qza"
FILTERED_REPSEQS="${QIIME_CORE}/rep-seqs_min${MIN_SAMPLES_FEATURE}samples.qza"
[[ -f "${FILTERED_TABLE}" ]] || FILTERED_TABLE="${QIIME_CORE}/table.qza"
[[ -f "${FILTERED_REPSEQS}" ]] || FILTERED_REPSEQS="${QIIME_CORE}/rep-seqs.qza"

# 6. Arbre : MAFFT -> mask -> FastTree -> enracinement midpoint.
if [[ "${RUN_TREE}" == true ]]; then
    log "Construction de l'arbre phylogenetique des ASV filtrees"
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

# 7. Exports utiles (BIOM, FASTA, TSV de statistiques, Newick).
log "Export des resultats QIIME 2"
mkdir -p "${QIIME_EXPORT}/core" "${QIIME_EXPORT}/tree" "${QIIME_EXPORT}/visual"
qiime tools export --input-path "${QIIME_CORE}/table.qza" --output-path "${QIIME_EXPORT}/core/table"
qiime tools export --input-path "${QIIME_CORE}/rep-seqs.qza" --output-path "${QIIME_EXPORT}/core/rep-seqs"
qiime tools export --input-path "${QIIME_CORE}/denoising-stats.qza" --output-path "${QIIME_EXPORT}/core/denoising-stats"
if [[ -f "${QIIME_CORE}/table_min${MIN_SAMPLES_FEATURE}samples.qza" ]]; then
    qiime tools export --input-path "${QIIME_CORE}/table_min${MIN_SAMPLES_FEATURE}samples.qza" --output-path "${QIIME_EXPORT}/core/table_min${MIN_SAMPLES_FEATURE}samples"
    qiime tools export --input-path "${QIIME_CORE}/rep-seqs_min${MIN_SAMPLES_FEATURE}samples.qza" --output-path "${QIIME_EXPORT}/core/rep-seqs_min${MIN_SAMPLES_FEATURE}samples"
fi
if [[ -f "${QIIME_TREE}/rooted-tree.qza" ]]; then
    qiime tools export --input-path "${QIIME_TREE}/aligned-rep-seqs.qza" --output-path "${QIIME_EXPORT}/tree/aligned-rep-seqs"
    qiime tools export --input-path "${QIIME_TREE}/masked-aligned-rep-seqs.qza" --output-path "${QIIME_EXPORT}/tree/masked-aligned-rep-seqs"
    qiime tools export --input-path "${QIIME_TREE}/unrooted-tree.qza" --output-path "${QIIME_EXPORT}/tree/unrooted-tree"
    qiime tools export --input-path "${QIIME_TREE}/rooted-tree.qza" --output-path "${QIIME_EXPORT}/tree/rooted-tree"
fi

for qzv in "${QIIME_VISUAL}"/*.qzv; do
    [[ -e "${qzv}" ]] || continue
    name="$(basename "${qzv}" .qzv)"
    qiime tools export --input-path "${qzv}" --output-path "${QIIME_EXPORT}/visual/${name}"
done

log "Pipeline termine avec succes."
log "Manifest : ${MANIFEST}"
log "Metadata : ${METADATA}"
log "Resultats QIIME 2 : ${QIIME_DIR}"
