#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

############################################
# JEDI cross-domain pipeline for Valormicro
# Complete fixed version v5
############################################

SCRIPT_VERSION="JEDI_VALORMICRO_FIXED_V5_2026-10-02"

PROJECT_ROOT="${PROJECT_ROOT:-/nvme/bio/data_fungi/valormicro_V2}"
RAW_DIR="${PROJECT_ROOT}/01_raw_data"
OLD_PIPE_DIR="${PROJECT_ROOT}/02_amplicon_pipeline"
OLD_SAMPLES_TSV="${OLD_PIPE_DIR}/04_database_files/samples.tsv"
XLSX_METADATA="${RAW_DIR}/00_infos_data.xlsx"

JEDI_DIR="${PROJECT_ROOT}/03_JEDI_pipeline"
INPUT_DIR="${JEDI_DIR}/00_inputs"
SILVA_DIR="${JEDI_DIR}/01_nfcore_silva"
PR2_DIR="${JEDI_DIR}/02_nfcore_pr2"
INTEGRATED_DIR="${JEDI_DIR}/03_integrated"
REF_DIR="${JEDI_DIR}/reference_databases"
LOG_DIR="${JEDI_DIR}/logs"
TMP_DIR="${JEDI_DIR}/tmp"

CONTAINER_ROOT="${JEDI_DIR}/container_cache"
SINGULARITY_TMPDIR_LOCAL="${CONTAINER_ROOT}/singularity_tmp"
SINGULARITY_LAYER_CACHE="${CONTAINER_ROOT}/singularity_layers"
SINGULARITY_IMAGE_CACHE="${CONTAINER_ROOT}/singularity_images"

LOCK_DIR="${JEDI_DIR}/.jedi_pipeline.lock"
SUCCESS_FLAG="${INTEGRATED_DIR}/PIPELINE_SUCCESS.txt"
FAIL_FLAG="${INTEGRATED_DIR}/PIPELINE_FAILED.txt"
RUN_MANIFEST="${INTEGRATED_DIR}/run_manifest.tsv"
LOG_FILE="${LOG_DIR}/jedi_pipeline_fixed_v5.log"

SAMPLESHEET="${INPUT_DIR}/samplesheet_jedi.tsv"
SAMPLESHEET_DEBUG="${INPUT_DIR}/samplesheet_jedi_debug.tsv"
SILVA_PARAMS="${INPUT_DIR}/params_silva_fixed_v5.yaml"
NF_INFRA_CFG="${INPUT_DIR}/nextflow_infrastructure.config"

SILVA_ASV_FASTA="${SILVA_DIR}/dada2/ASV_seqs.fasta"
SILVA_ASV_TABLE="${SILVA_DIR}/dada2/ASV_table.tsv"
SILVA_TAXONOMY="${SILVA_DIR}/dada2/ASV_tax.silva_138_2.tsv"

PR2_R_SCRIPT="${PR2_DIR}/classify_pr2_from_silva_asvs.R"
PR2_TAXONOMY_TSV="${PR2_DIR}/taxonomy_PR2_ASV.tsv"
PR2_REFERENCE_FASTA="${REF_DIR}/pr2_version_5.1.0_SSU_dada2.fasta.gz"
PR2_REFERENCE_URL="https://github.com/pr2database/pr2database/releases/download/v5.1.0.0/pr2_version_5.1.0_SSU_dada2.fasta.gz"

NFCORE_VERSION="${NFCORE_VERSION:-2.18.0}"
PROFILE="${PROFILE:-singularity}"
AMPLISEQ_REPO="${AMPLISEQ_REPO:-nf-core/ampliseq}"

PYTHON_BIN="${PYTHON_BIN:-python3}"
R_BIN="${R_BIN:-Rscript}"
PR2_CONDA_ENV_CANDIDATES="${PR2_CONDA_ENV_CANDIDATES:-qiime2-amplicon-2024.10 dada2 qiime2-amplicon-2024.5}"

FORWARD_PRIMER="${FORWARD_PRIMER:-GTGYCAGCMGCCGCGGTAA}"
REVERSE_PRIMER="${REVERSE_PRIMER:-CCGYCAATTYMTTTRAGTTT}"
TRUNC_LEN_F="${TRUNC_LEN_F:-231}"
TRUNC_LEN_R="${TRUNC_LEN_R:-230}"
MAX_EE="${MAX_EE:-2}"
TRUNC_QMIN="${TRUNC_QMIN:-25}"
TRUNC_RMIN="${TRUNC_RMIN:-0.75}"
SILVA_REF="${SILVA_REF:-silva=138.2}"
PR2_REF_LABEL="PR2_v5.1.0_dada2"

PRELOAD_IMAGE_DOCKER="docker://biocontainers/biocontainers:v1.2.0_cv1"
PRELOAD_IMAGE_SIF="${SINGULARITY_IMAGE_CACHE}/biocontainers_v1.2.0_cv1.sif"
PRELOAD_IMAGE_ALIAS1="${SINGULARITY_IMAGE_CACHE}/containers.biocontainers.pro-s3-SingImgsRepo-biocontainers-v1.2.0_cv1-biocontainers_v1.2.0_cv1.img.img"
PRELOAD_IMAGE_ALIAS2="${CONTAINER_ROOT}/containers.biocontainers.pro-s3-SingImgsRepo-biocontainers-v1.2.0_cv1-biocontainers_v1.2.0_cv1.img.img"
PULL_TIMEOUT="${PULL_TIMEOUT:-12 h}"

FORCE_SILVA="${FORCE_SILVA:-0}"
FORCE_PR2="${FORCE_PR2:-0}"
FORCE_INTEGRATION="${FORCE_INTEGRATION:-0}"
FORCE_ALL="${FORCE_ALL:-0}"
KEEP_FAILURE_MARKER="${KEEP_FAILURE_MARKER:-0}"

mkdir -p \
  "${INPUT_DIR}" "${SILVA_DIR}" "${PR2_DIR}" "${INTEGRATED_DIR}" \
  "${REF_DIR}" "${LOG_DIR}" "${TMP_DIR}" "${CONTAINER_ROOT}" \
  "${SINGULARITY_TMPDIR_LOCAL}" "${SINGULARITY_LAYER_CACHE}" \
  "${SINGULARITY_IMAGE_CACHE}" "${INTEGRATED_DIR}/tables_by_domain"

touch "${LOG_FILE}"
exec > >(tee -a "${LOG_FILE}") 2>&1

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
  log "ERREUR : $*"
  printf 'FAILED\t%s\t%s\n' "$(date '+%F %T')" "$*" > "${FAIL_FLAG}" || true
  exit 1
}

on_error() {
  local exit_code="$?"
  local line_no="${1:-unknown}"
  local cmd="${2:-unknown}"
  trap - ERR
  log "ERREUR : code ${exit_code}, ligne ${line_no}, commande : ${cmd}"
  printf 'FAILED\t%s\tline=%s\tcmd=%s\texit=%s\n' \
    "$(date '+%F %T')" "${line_no}" "${cmd}" "${exit_code}" > "${FAIL_FLAG}" || true
  exit "${exit_code}"
}
trap 'on_error ${LINENO} "${BASH_COMMAND}"' ERR

cleanup_lock() {
  rm -rf "${LOCK_DIR}" || true
}
trap cleanup_lock EXIT

acquire_lock() {
  if mkdir "${LOCK_DIR}" 2>/dev/null; then
    printf '%s\n' "$$" > "${LOCK_DIR}/pid"
    return 0
  fi

  local old_pid=""
  if [[ -s "${LOCK_DIR}/pid" ]]; then
    old_pid="$(cat "${LOCK_DIR}/pid" 2>/dev/null || true)"
  fi
  if [[ -n "${old_pid}" ]] && kill -0 "${old_pid}" 2>/dev/null; then
    die "Un autre run JEDI est actif (PID ${old_pid}) : ${LOCK_DIR}"
  fi

  log "Suppression d'un verrou obsolète : ${LOCK_DIR}"
  rm -rf "${LOCK_DIR}"
  mkdir "${LOCK_DIR}"
  printf '%s\n' "$$" > "${LOCK_DIR}/pid"
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Commande introuvable : $1"
}

setup_environment() {
  export NXF_SINGULARITY_CACHEDIR="${SINGULARITY_IMAGE_CACHE}"
  export NXF_APPTAINER_CACHEDIR="${SINGULARITY_IMAGE_CACHE}"
  export SINGULARITY_CACHEDIR="${SINGULARITY_LAYER_CACHE}"
  export APPTAINER_CACHEDIR="${SINGULARITY_LAYER_CACHE}"
  export SINGULARITY_TMPDIR="${SINGULARITY_TMPDIR_LOCAL}"
  export APPTAINER_TMPDIR="${SINGULARITY_TMPDIR_LOCAL}"
  export TMPDIR="${SINGULARITY_TMPDIR_LOCAL}"
  export TEMP="${SINGULARITY_TMPDIR_LOCAL}"
  export TMP="${SINGULARITY_TMPDIR_LOCAL}"
  export NXF_HOME="${JEDI_DIR}/.nextflow"
  export NXF_OFFLINE='false'
  mkdir -p "${NXF_HOME}" "${SINGULARITY_TMPDIR_LOCAL}" \
    "${SINGULARITY_LAYER_CACHE}" "${SINGULARITY_IMAGE_CACHE}"
}

activate_nextflow_java() {
  local conda_base=""
  if command -v conda >/dev/null 2>&1; then
    conda_base="$(conda info --base 2>/dev/null || true)"
    if [[ -n "${conda_base}" && -f "${conda_base}/etc/profile.d/conda.sh" ]]; then
      # shellcheck disable=SC1090
      source "${conda_base}/etc/profile.d/conda.sh"
      conda activate nextflow-java || true
    fi
  fi

  require_cmd java
  require_cmd nextflow
  export JAVA_HOME="$(dirname "$(dirname "$(readlink -f "$(command -v java)")")")"
  export NXF_JAVA_HOME="${JAVA_HOME}"
  unset JAVA_CMD || true
  hash -r
  java -version
  nextflow -version
}

check_prerequisites() {
  require_cmd bash
  require_cmd awk
  require_cmd sed
  require_cmd grep
  require_cmd gzip
  require_cmd find
  require_cmd singularity
  require_cmd nextflow
  require_cmd curl
  require_cmd "${PYTHON_BIN}"
  require_cmd java

  [[ -d "${RAW_DIR}" ]] || die "Répertoire FASTQ absent : ${RAW_DIR}"
  [[ -s "${OLD_SAMPLES_TSV}" || -s "${XLSX_METADATA}" ]] || die \
    "Aucune source de métadonnées : ${OLD_SAMPLES_TSV} ou ${XLSX_METADATA}"
}

write_nextflow_config() {
  cat > "${NF_INFRA_CFG}" <<EOF_CFG
process {
  executor = 'local'
  scratch = false
}

singularity {
  enabled = true
  autoMounts = true
  cacheDir = '${SINGULARITY_IMAGE_CACHE}'
  pullTimeout = '${PULL_TIMEOUT}'
}

apptainer.enabled = false
docker.enabled = false
cleanup = false
report.enabled = false
timeline.enabled = false
trace.enabled = true
dag.enabled = false
EOF_CFG
}

preload_problematic_container() {
  log "Préchargement sécurisé du conteneur BioContainers problématique"
  rm -f "${PRELOAD_IMAGE_ALIAS1}.pulling."* "${PRELOAD_IMAGE_ALIAS2}.pulling."* 2>/dev/null || true
  if [[ ! -s "${PRELOAD_IMAGE_SIF}" ]]; then
    log "Téléchargement local du conteneur depuis ${PRELOAD_IMAGE_DOCKER}"
    singularity pull "${PRELOAD_IMAGE_SIF}" "${PRELOAD_IMAGE_DOCKER}"
  else
    log "Conteneur déjà présent : ${PRELOAD_IMAGE_SIF}"
  fi
  singularity exec "${PRELOAD_IMAGE_SIF}" python --version >/dev/null
  ln -sfn "${PRELOAD_IMAGE_SIF}" "${PRELOAD_IMAGE_ALIAS1}"
  ln -sfn "${PRELOAD_IMAGE_SIF}" "${PRELOAD_IMAGE_ALIAS2}"
  log "Alias Nextflow créés vers ${PRELOAD_IMAGE_SIF}"
}

build_samplesheet() {
  log "Construction et validation du samplesheet : ${SAMPLESHEET}"

  "${PYTHON_BIN}" - \
    "${OLD_SAMPLES_TSV}" "${RAW_DIR}" "${SAMPLESHEET}" \
    "${SAMPLESHEET_DEBUG}" "${XLSX_METADATA}" <<'PY_SAMPLESHEET'
import csv
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

old_tsv, raw_dir, out_tsv, debug_tsv, xlsx = sys.argv[1:6]


def nfcore_id(value):
    value = str(value).strip()
    value = re.sub(r'[^A-Za-z0-9_]+', '_', value)
    value = re.sub(r'_+', '_', value).strip('_')
    if not value:
        value = 'Sample'
    if not re.match(r'^[A-Za-z]', value):
        value = 'S_' + value
    if not re.match(r'^[A-Za-z][A-Za-z0-9_]*$', value):
        raise SystemExit('ID nf-core invalide : {}'.format(value))
    return value


def resolve_fastq(root_dir, filename):
    filename = os.path.basename(str(filename).strip())
    if not filename:
        return None
    direct = os.path.join(root_dir, filename)
    if os.path.isfile(direct):
        return os.path.abspath(direct)
    found = []
    for root, _, files in os.walk(root_dir):
        if filename in files:
            found.append(os.path.abspath(os.path.join(root, filename)))
    if len(found) == 1:
        return found[0]
    if not found:
        return None
    raise SystemExit('FASTQ ambigu pour {}:\n{}'.format(filename, '\n'.join(found)))


def write_outputs(records, debug_records):
    with open(out_tsv, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample', 'fastq_1', 'fastq_2'])
        for rec in records:
            writer.writerow([rec['sample'], rec['fastq_1'], rec['fastq_2']])
    with open(debug_tsv, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample', 'fastq_1', 'fastq_2', 'source'])
        for rec in debug_records:
            writer.writerow([rec['sample'], rec['fastq_1'], rec['fastq_2'], rec['source']])


def build_from_old_tsv(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None
    records, debug_records, used = [], [], set()
    with open(path, 'r', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames:
            raise SystemExit('En-tête absent dans {}'.format(path))
        fields = {name.strip().lower(): name for name in reader.fieldnames if name}

        def field(*candidates):
            for candidate in candidates:
                if candidate.lower() in fields:
                    return fields[candidate.lower()]
            return None

        r1_col = field('R1', 'fastq_1', 'fastq1')
        r2_col = field('R2', 'fastq_2', 'fastq2')
        label_col = field('New_label', 'Newlabel', 'sample-id', 'sample', 'Label')
        if not all((r1_col, r2_col, label_col)):
            raise SystemExit('Colonnes R1/R2/label absentes dans {}'.format(path))

        for line_num, row in enumerate(reader, start=2):
            r1_name = (row.get(r1_col) or '').strip()
            r2_name = (row.get(r2_col) or '').strip()
            label = (row.get(label_col) or '').strip()
            if not r1_name and not r2_name and not label:
                continue
            if not all((r1_name, r2_name, label)):
                raise SystemExit('Ligne {} incomplète dans {}'.format(line_num, path))
            sample = nfcore_id(label)
            if sample in used:
                raise SystemExit('ID dupliqué après normalisation : {}'.format(sample))
            r1_abs = resolve_fastq(raw_dir, r1_name)
            r2_abs = resolve_fastq(raw_dir, r2_name)
            if not r1_abs or not r2_abs:
                raise SystemExit("FASTQ introuvable ligne {} : R1='{}', R2='{}'".format(
                    line_num, r1_name, r2_name))
            used.add(sample)
            rec = {'sample': sample, 'fastq_1': r1_abs, 'fastq_2': r2_abs}
            records.append(rec)
            debug_records.append(dict(rec, source='samples.tsv'))
    return records, debug_records


def header_key(value):
    value = str(value).replace('\xa0', ' ').strip().lower()
    return re.sub(r'[^a-z0-9]', '', value)


def read_xlsx(path):
    ns = {'a': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}

    def col_to_idx(col):
        n = 0
        for char in col:
            if char.isalpha():
                n = n * 26 + ord(char.upper()) - 64
        return n - 1

    with zipfile.ZipFile(path) as archive:
        strings = []
        if 'xl/sharedStrings.xml' in archive.namelist():
            root = ET.fromstring(archive.read('xl/sharedStrings.xml'))
            for si in root.findall('a:si', ns):
                strings.append(''.join((x.text or '') for x in si.iterfind('.//a:t', ns)))

        workbook = ET.fromstring(archive.read('xl/workbook.xml'))
        rel_root = ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))
        rels = {rel.attrib['Id']: rel.attrib['Target'] for rel in rel_root}
        sheets = []
        for sheet in workbook.find('a:sheets', ns):
            rid = sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
            target = rels[rid].lstrip('/')
            if not target.startswith('xl/'):
                target = 'xl/' + target
            sheets.append((sheet.attrib['name'], target))

        all_rows = []
        for sheet_name, target in sheets:
            root = ET.fromstring(archive.read(target))
            sheet_data = root.find('a:sheetData', ns)
            rows = []
            if sheet_data is None:
                all_rows.append((sheet_name, rows))
                continue
            for row in sheet_data.findall('a:row', ns):
                values = {}
                for cell in row.findall('a:c', ns):
                    ref = cell.attrib.get('r', 'A1')
                    idx = col_to_idx(''.join(c for c in ref if c.isalpha()))
                    cell_type = cell.attrib.get('t')
                    value_node = cell.find('a:v', ns)
                    value = ''
                    if value_node is not None and value_node.text is not None:
                        value = value_node.text
                        if cell_type == 's':
                            value = strings[int(value)]
                    else:
                        inline = cell.find('a:is', ns)
                        if inline is not None:
                            value = ''.join((x.text or '') for x in inline.iterfind('.//a:t', ns))
                    values[idx] = str(value).strip()
                if values:
                    row_values = [''] * (max(values) + 1)
                    for idx, value in values.items():
                        row_values[idx] = value
                    rows.append(row_values)
            all_rows.append((sheet_name, rows))
        return all_rows


def build_from_xlsx(path):
    if not os.path.isfile(path):
        raise SystemExit('Fichier Excel introuvable : {}'.format(path))
    records, debug_records, used = [], [], set()
    for sheet_name, rows in read_xlsx(path):
        header = None
        header_idx = None
        for idx, row in enumerate(rows[:20]):
            keys = [header_key(x) for x in row]
            if 'r1' in keys and 'r2' in keys:
                header, header_idx = row, idx
                break
        if header is None:
            continue

        keys = [header_key(x) for x in header]

        def find_col(names):
            accepted = {header_key(x) for x in names}
            for idx, key in enumerate(keys):
                if key in accepted:
                    return idx
            return None

        r1_col = find_col(('R1', 'Read1', 'Forward', 'Forward read'))
        r2_col = find_col(('R2', 'Read2', 'Reverse', 'Reverse read'))
        label_col = find_col(('Newlabel', 'New label', 'New_label', 'Sample',
                              'SampleID', 'Sample ID', 'Echantillon', 'Échantillon', 'Label'))
        if r1_col is None or r2_col is None or label_col is None:
            continue

        for row_num, row in enumerate(rows[header_idx + 1:], start=header_idx + 2):
            if not any(str(x).strip() for x in row):
                continue
            r1_name = row[r1_col].strip() if r1_col < len(row) else ''
            r2_name = row[r2_col].strip() if r2_col < len(row) else ''
            label = row[label_col].strip() if label_col < len(row) else ''
            if not r1_name and not r2_name and not label:
                continue
            if not all((r1_name, r2_name, label)):
                raise SystemExit("Feuille '{}', ligne {} incomplète".format(sheet_name, row_num))
            sample = nfcore_id(label)
            if sample in used:
                raise SystemExit('ID dupliqué après normalisation : {}'.format(sample))
            r1_abs = resolve_fastq(raw_dir, r1_name)
            r2_abs = resolve_fastq(raw_dir, r2_name)
            if not r1_abs or not r2_abs:
                raise SystemExit("Feuille '{}', ligne {} : FASTQ introuvable(s)".format(
                    sheet_name, row_num))
            used.add(sample)
            rec = {'sample': sample, 'fastq_1': r1_abs, 'fastq_2': r2_abs}
            records.append(rec)
            debug_records.append(dict(rec, source='00_infos_data.xlsx::{}'.format(sheet_name)))

    if not records:
        raise SystemExit("Aucune paire valide créée depuis {}".format(path))
    return records, debug_records


result = build_from_old_tsv(old_tsv)
if result is None:
    result = build_from_xlsx(xlsx)
records, debug_records = result
if not records:
    raise SystemExit('Aucun échantillon valide détecté')
write_outputs(records, debug_records)
print('OK : {} échantillons valides écrits.'.format(len(records)))
PY_SAMPLESHEET

  sed -i 's/\r$//' "${SAMPLESHEET}" "${SAMPLESHEET_DEBUG}"
  local sample_count
  sample_count="$(awk 'END {print NR - 1}' "${SAMPLESHEET}")"
  [[ "${sample_count}" -gt 0 ]] || die "Samplesheet vide"
  log "${sample_count} échantillons écrits dans ${SAMPLESHEET}"
}

check_fastq_integrity() {
  log "Contrôle gzip de tous les FASTQ"
  local count=0
  local sample r1 r2
  while IFS=$'\t' read -r sample r1 r2; do
    [[ "${sample}" == "sample" ]] && continue
    sample="${sample//$'\r'/}"
    r1="${r1//$'\r'/}"
    r2="${r2//$'\r'/}"
    [[ -f "${r1}" ]] || die "FASTQ R1 absent : ${r1}"
    [[ -f "${r2}" ]] || die "FASTQ R2 absent : ${r2}"
    gzip -t "${r1}"
    gzip -t "${r2}"
    count=$((count + 1))
  done < <(tr -d '\r' < "${SAMPLESHEET}")
  log "${count} paires FASTQ validées"
}

write_silva_params() {
  cat > "${SILVA_PARAMS}" <<EOF_SILVA
input: "${SAMPLESHEET}"
outdir: "${SILVA_DIR}"

FW_primer: "${FORWARD_PRIMER}"
RV_primer: "${REVERSE_PRIMER}"
retain_untrimmed: false
cutadapt_min_overlap: 5
cutadapt_max_error_rate: 0.1

truncq: 2
trunclenf: ${TRUNC_LEN_F}
trunclenr: ${TRUNC_LEN_R}
trunc_qmin: ${TRUNC_QMIN}
trunc_rmin: ${TRUNC_RMIN}
max_ee: ${MAX_EE}
min_len: 50

single_end: false
multiple_sequencing_runs: false
sample_inference: "independent"
mergepairs_strategy: "consensus"
mergepairs_consensus_match: 1
mergepairs_consensus_mismatch: -2
mergepairs_consensus_gap: -4
mergepairs_consensus_minoverlap: 12
mergepairs_consensus_maxmismatch: 0
mergepairs_consensus_percentile_cutoff: 0.001

dada_ref_taxonomy: "${SILVA_REF}"
cut_dada_ref_taxonomy: true
exclude_taxa: "none"

skip_fastqc: false
skip_cutadapt: false
skip_barrnap: false
skip_qiime: true
skip_dada_addspecies: true
skip_abundance_tables: true
skip_alpha_rarefaction: true
skip_diversity_indices: true

save_intermediates: true
report_title: "JEDI valormicro - ASV et SILVA"
EOF_SILVA
}

run_nfcore() {
  local params_file="$1"
  local work_subdir="$2"
  local label="$3"
  log "Lancement ${AMPLISEQ_REPO} ${label}"
  log "Paramètres : ${params_file}"
  log "Work : ${work_subdir}"
  (
    cd "${JEDI_DIR}"
    nextflow run "${AMPLISEQ_REPO}" \
      -r "${NFCORE_VERSION}" \
      -profile "${PROFILE}" \
      -params-file "${params_file}" \
      -work-dir "${work_subdir}" \
      -c "${NF_INFRA_CFG}" \
      -resume \
      -ansi-log false
  )
}

has_silva_success() {
  [[ -s "${SILVA_ASV_FASTA}" ]] \
    && [[ -s "${SILVA_ASV_TABLE}" ]] \
    && [[ -s "${SILVA_TAXONOMY}" ]]
}

run_silva_pipeline() {
  if [[ "${FORCE_ALL}" == "1" || "${FORCE_SILVA}" == "1" ]]; then
    log "Forçage du recalcul SILVA"
    rm -rf "${SILVA_DIR:?}/"*
  fi

  if has_silva_success; then
    log "ÉTAPE 1/3 SILVA déjà disponible, réutilisation"
  else
    log "ÉTAPE 1/3 DADA2/JEDI et taxonomie SILVA"
    write_silva_params
    run_nfcore "${SILVA_PARAMS}" "${JEDI_DIR}/work_silva" "SILVA"
  fi

  [[ -s "${SILVA_ASV_FASTA}" ]] || die "FASTA ASV SILVA absent : ${SILVA_ASV_FASTA}"
  [[ -s "${SILVA_ASV_TABLE}" ]] || die "Table ASV SILVA absente : ${SILVA_ASV_TABLE}"
  [[ -s "${SILVA_TAXONOMY}" ]] || die "Taxonomie SILVA absente : ${SILVA_TAXONOMY}"
}

download_pr2_reference() {
  if [[ -s "${PR2_REFERENCE_FASTA}" ]]; then
    gzip -t "${PR2_REFERENCE_FASTA}" || die "Référence PR2 gzip corrompue : ${PR2_REFERENCE_FASTA}"
    log "Référence PR2 déjà présente : ${PR2_REFERENCE_FASTA}"
    return 0
  fi

  log "Téléchargement de la référence PR2 DADA2"
  local tmp_ref="${PR2_REFERENCE_FASTA}.part"
  rm -f "${tmp_ref}"
  curl --fail --location --retry 3 --retry-delay 5 \
    --output "${tmp_ref}" "${PR2_REFERENCE_URL}"
  gzip -t "${tmp_ref}" || die "Référence PR2 téléchargée mais corrompue"
  mv "${tmp_ref}" "${PR2_REFERENCE_FASTA}"
}

write_pr2_r_script() {
  cat > "${PR2_R_SCRIPT}" <<'RSCRIPT'
suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: script.R ASV.fasta PR2.fasta.gz output.tsv")
}

asv_fasta <- args[[1]]
pr2_fasta <- args[[2]]
out_tsv <- args[[3]]

read_fasta_named <- function(path) {
  x <- readLines(path, warn = FALSE)
  headers <- which(startsWith(x, ">"))
  if (length(headers) == 0L) {
    stop("Aucune entrée FASTA dans : ", path)
  }
  ends <- c(headers[-1L] - 1L, length(x))
  ids <- sub("^>", "", x[headers])
  ids <- sub("[[:space:]].*$", "", ids)
  seqs <- vapply(seq_along(headers), function(i) {
    paste0(x[(headers[i] + 1L):ends[i]], collapse = "")
  }, character(1))
  seqs <- toupper(seqs)
  if (anyDuplicated(ids)) {
    stop("Identifiants FASTA dupliqués")
  }
  names(seqs) <- ids
  seqs
}

seqs <- read_fasta_named(asv_fasta)
message("Classification PR2 de ", length(seqs), " ASV")

tax <- assignTaxonomy(
  seqs = unname(seqs),
  refFasta = pr2_fasta,
  minBoot = 50,
  tryRC = TRUE,
  outputBootstraps = FALSE,
  multithread = TRUE
)

tax <- as.data.frame(tax, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(tax) != length(seqs)) {
  stop("Nombre de classifications différent du nombre d'ASV")
}
tax$ASV <- names(seqs)

rename_map <- c(
  Kingdom = "domain",
  Domain = "domain",
  Supergroup = "supergroup",
  Division = "division",
  Phylum = "division",
  Class = "class",
  Order = "order",
  Family = "family",
  Genus = "genus",
  Species = "species"
)
for (old in names(rename_map)) {
  if (old %in% colnames(tax)) {
    colnames(tax)[colnames(tax) == old] <- rename_map[[old]]
  }
}

rank_order <- c(
  "ASV", "domain", "supergroup", "division",
  "class", "order", "family", "genus", "species"
)
for (rank in setdiff(rank_order, colnames(tax))) {
  tax[[rank]] <- NA_character_
}
tax <- tax[, rank_order, drop = FALSE]

write.table(
  tax,
  file = out_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)
RSCRIPT
}

find_dada2_container() {
  local work_root="${JEDI_DIR}/work_silva"
  local cmd_file run_file image

  while IFS= read -r -d '' cmd_file; do
    grep -q 'library(dada2)' "${cmd_file}" 2>/dev/null || continue
    run_file="$(dirname "${cmd_file}")/.command.run"
    [[ -s "${run_file}" ]] || continue
    image="$(grep -oE '/[^[:space:]]+\.(sif|img)' "${run_file}" | tr -d '"' | tail -n 1 || true)"
    if [[ -n "${image}" && -s "${image}" ]]; then
      printf '%s\n' "${image}"
      return 0
    fi
  done < <(find "${work_root}" -type f -name '.command.sh' -print0 2>/dev/null)

  return 1
}

run_pr2_classifier() {
  local conda_base=""
  local env_name=""
  local dada2_image=""

  if command -v conda >/dev/null 2>&1; then
    conda_base="$(conda info --base 2>/dev/null || true)"
    if [[ -n "${conda_base}" && -f "${conda_base}/etc/profile.d/conda.sh" ]]; then
      # shellcheck disable=SC1090
      source "${conda_base}/etc/profile.d/conda.sh"
      for env_name in ${PR2_CONDA_ENV_CANDIDATES}; do
        if conda activate "${env_name}" >/dev/null 2>&1 \
          && command -v "${R_BIN}" >/dev/null 2>&1 \
          && "${R_BIN}" --vanilla -e 'suppressPackageStartupMessages(library(dada2))' >/dev/null 2>&1; then
          log "Environnement R/dada2 activé pour PR2 : ${env_name}"
          "${R_BIN}" --vanilla "${PR2_R_SCRIPT}" \
            "${SILVA_ASV_FASTA}" "${PR2_REFERENCE_FASTA}" "${PR2_TAXONOMY_TSV}"
          return 0
        fi
      done
    fi
  fi

  if command -v "${R_BIN}" >/dev/null 2>&1 \
    && "${R_BIN}" --vanilla -e 'suppressPackageStartupMessages(library(dada2))' >/dev/null 2>&1; then
    log "Utilisation du R/dada2 courant pour PR2"
    "${R_BIN}" --vanilla "${PR2_R_SCRIPT}" \
      "${SILVA_ASV_FASTA}" "${PR2_REFERENCE_FASTA}" "${PR2_TAXONOMY_TSV}"
    return 0
  fi

  dada2_image="$(find_dada2_container || true)"
  if [[ -n "${dada2_image}" ]]; then
    log "Utilisation du conteneur DADA2 nf-core pour PR2 : ${dada2_image}"
    singularity exec --bind "${PROJECT_ROOT}:${PROJECT_ROOT}" "${dada2_image}" \
      Rscript --vanilla "${PR2_R_SCRIPT}" \
      "${SILVA_ASV_FASTA}" "${PR2_REFERENCE_FASTA}" "${PR2_TAXONOMY_TSV}"
    return 0
  fi

  die "Aucun environnement R/dada2 utilisable pour la classification PR2"
}

classify_pr2_from_silva_asvs() {
  [[ -s "${SILVA_ASV_FASTA}" ]] || die \
    "FASTA ASV SILVA absent pour PR2 : ${SILVA_ASV_FASTA}"

  if [[ "${FORCE_ALL}" == "1" || "${FORCE_PR2}" == "1" ]]; then
    log "Forçage du recalcul PR2"
    rm -f "${PR2_TAXONOMY_TSV}" "${PR2_R_SCRIPT}"
  fi

  if [[ -s "${PR2_TAXONOMY_TSV}" ]]; then
    log "ÉTAPE 2/3 PR2 déjà disponible, réutilisation"
    return 0
  fi

  mkdir -p "${PR2_DIR}"
  download_pr2_reference
  write_pr2_r_script
  log "ÉTAPE 2/3 classification PR2 directe des ASV SILVA"
  run_pr2_classifier
  [[ -s "${PR2_TAXONOMY_TSV}" ]] || die \
    "Taxonomie PR2 non produite : ${PR2_TAXONOMY_TSV}"
}

integrate_outputs() {
  log "ÉTAPE 3/3 intégration SILVA + PR2"

  if [[ "${FORCE_ALL}" == "1" || "${FORCE_INTEGRATION}" == "1" ]]; then
    log "Forçage de l'intégration finale"
    rm -f "${INTEGRATED_DIR}"/*.tsv "${INTEGRATED_DIR}/tables_by_domain"/*.tsv 2>/dev/null || true
  fi

  [[ -s "${SILVA_ASV_FASTA}" ]] || die "FASTA ASV absent : ${SILVA_ASV_FASTA}"
  [[ -s "${SILVA_ASV_TABLE}" ]] || die "Table ASV absente : ${SILVA_ASV_TABLE}"
  [[ -s "${SILVA_TAXONOMY}" ]] || die "Taxonomie SILVA absente : ${SILVA_TAXONOMY}"
  [[ -s "${PR2_TAXONOMY_TSV}" ]] || die "Taxonomie PR2 absente : ${PR2_TAXONOMY_TSV}"

  mkdir -p "${INTEGRATED_DIR}/tables_by_domain"

  SILVA_TABLE="${SILVA_ASV_TABLE}" \
  SILVA_ASV="${SILVA_ASV_FASTA}" \
  SILVA_TAX="${SILVA_TAXONOMY}" \
  PR2_TAX="${PR2_TAXONOMY_TSV}" \
  INTEGRATED_DIR="${INTEGRATED_DIR}" \
  "${PYTHON_BIN}" <<'PY_INTEGRATE'
import csv
import math
import os
import re
from collections import defaultdict, OrderedDict

integrated_dir = os.environ['INTEGRATED_DIR']
silva_table = os.environ['SILVA_TABLE']
silva_asv = os.environ['SILVA_ASV']
silva_tax_path = os.environ['SILVA_TAX']
pr2_tax_path = os.environ['PR2_TAX']

os.makedirs(integrated_dir, exist_ok=True)
os.makedirs(os.path.join(integrated_dir, 'tables_by_domain'), exist_ok=True)

RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def parse_fasta(path):
    sequences = OrderedDict()
    current, buffer = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current is not None:
                    sequences[current] = ''.join(buffer)
                current = line[1:].split()[0]
                buffer = []
            else:
                buffer.append(line)
    if current is not None:
        sequences[current] = ''.join(buffer)
    if not sequences:
        raise SystemExit('Aucune séquence dans {}'.format(path))
    return sequences


def sniff_delimiter(line):
    return '\t' if line.count('\t') >= line.count(',') else ','


def normalize_header(value):
    value = (value or '').strip().lower().replace('kingdom', 'domain')
    return re.sub(r'[^a-z0-9]+', '_', value).strip('_')


def clean_tax(value):
    value = (value or '').strip()
    if value.lower() in ('', 'na', 'nan', 'none', 'unassigned', 'unknown'):
        return ''
    return value


def parse_taxonomy(path):
    with open(path) as handle:
        first = handle.readline()
    if not first:
        raise SystemExit('Table taxonomique vide : {}'.format(path))
    delimiter = sniff_delimiter(first)
    with open(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        rows = list(reader)
    if len(rows) < 2:
        raise SystemExit('Table taxonomique sans données : {}'.format(path))

    headers = [normalize_header(x) for x in rows[0]]
    id_aliases = {'asv', 'asv_id', 'featureid', 'feature_id', 'feature', 'otu', 'otu_id', 'id'}
    id_idx = next((i for i, h in enumerate(headers) if h in id_aliases), 0)
    aliases = {
        'domain': {'domain'},
        'phylum': {'phylum', 'division'},
        'class': {'class'},
        'order': {'order'},
        'family': {'family'},
        'genus': {'genus'},
        'species': {'species', 'species_exact'},
    }
    rank_idx = {}
    for rank, accepted in aliases.items():
        for idx, header in enumerate(headers):
            if header in accepted:
                rank_idx[rank] = idx
                break

    taxonomy_idx = next((i for i, h in enumerate(headers)
                         if h in ('tax', 'taxon', 'taxonomy')), None)
    result = {}
    for row in rows[1:]:
        if not row or id_idx >= len(row):
            continue
        asv = row[id_idx].strip()
        if not asv:
            continue
        ranks = {rank: '' for rank in RANKS}
        if rank_idx:
            for rank, idx in rank_idx.items():
                if idx < len(row):
                    ranks[rank] = clean_tax(row[idx])
        elif taxonomy_idx is not None and taxonomy_idx < len(row):
            parts = re.split(r'[;,]\s*', row[taxonomy_idx].strip())
            cleaned = [re.sub(r'^[dkpcofgs]__', '', part) for part in parts]
            for idx, rank in enumerate(RANKS):
                if idx < len(cleaned):
                    ranks[rank] = clean_tax(cleaned[idx])
        result[asv] = ranks
    return result


def parse_count_table(path):
    with open(path) as handle:
        first = handle.readline()
    delimiter = sniff_delimiter(first)
    counts = OrderedDict()
    with open(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        samples = [x.strip() for x in header[1:]]
        if not samples:
            raise SystemExit('Aucun échantillon dans {}'.format(path))
        for row in reader:
            if not row:
                continue
            asv = row[0].strip()
            if not asv:
                continue
            values = []
            for value in row[1:1 + len(samples)]:
                value = value.strip()
                values.append(int(float(value)) if value not in ('', 'NA', 'nan') else 0)
            if len(values) < len(samples):
                values.extend([0] * (len(samples) - len(values)))
            counts[asv] = dict(zip(samples, values))
    if not counts:
        raise SystemExit('Aucun ASV dans {}'.format(path))
    return samples, counts


def canonical_domain(value):
    value = clean_tax(value)
    lower = value.lower()
    if lower in ('bacteria', 'bacteriae'):
        return 'Bacteria'
    if lower in ('archaea', 'archaeaota'):
        return 'Archaea'
    if lower in ('eukaryota', 'eukaryote', 'eukarya'):
        return 'Eukaryota'
    return value


def choose_assignment(silva_row, pr2_row):
    silva_domain = canonical_domain(silva_row.get('domain', ''))
    pr2_domain = canonical_domain(pr2_row.get('domain', ''))
    silva_text = ' '.join(clean_tax(v) for v in silva_row.values()).lower()

    if 'mitochond' in silva_text:
        return 'Mitochondria', silva_row, 'SILVA'
    if 'chloroplast' in silva_text or 'plastid' in silva_text:
        return 'Plastid', silva_row, 'SILVA'
    if silva_domain in ('Bacteria', 'Archaea'):
        return silva_domain, silva_row, 'SILVA'
    if pr2_domain == 'Eukaryota':
        return 'Eukaryota', pr2_row, 'PR2'
    if silva_domain == 'Eukaryota':
        base = pr2_row if any(clean_tax(v) for v in pr2_row.values()) else silva_row
        source = 'PR2' if base is pr2_row else 'SILVA'
        return 'Eukaryota', base, source
    if silva_domain:
        return silva_domain, silva_row, 'SILVA'
    return 'Unresolved', ({rank: '' for rank in RANKS}), 'UNRESOLVED'


sequences = parse_fasta(silva_asv)
silva_tax = parse_taxonomy(silva_tax_path)
pr2_tax = parse_taxonomy(pr2_tax_path)
samples, counts = parse_count_table(silva_table)

count_ids = set(counts)
if not count_ids.intersection(sequences):
    raise SystemExit('Aucun identifiant commun entre ASV_table.tsv et ASV_seqs.fasta')
if not count_ids.intersection(silva_tax):
    raise SystemExit('Aucun identifiant commun entre ASV_table.tsv et la taxonomie SILVA')
if not count_ids.intersection(pr2_tax):
    raise SystemExit('Aucun identifiant commun entre ASV_table.tsv et la taxonomie PR2')

consensus = []
domain_summary = defaultdict(int)
domain_per_sample = {sample: defaultdict(int) for sample in samples}

with open(os.path.join(integrated_dir, 'taxonomy_SILVA_original.tsv'), 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['ASV'] + RANKS)
    for asv in counts:
        row = silva_tax.get(asv, {})
        writer.writerow([asv] + [clean_tax(row.get(rank, '')) for rank in RANKS])

with open(os.path.join(integrated_dir, 'taxonomy_PR2_original.tsv'), 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['ASV'] + RANKS)
    for asv in counts:
        row = pr2_tax.get(asv, {})
        writer.writerow([asv] + [clean_tax(row.get(rank, '')) for rank in RANKS])

with open(os.path.join(integrated_dir, 'taxonomy_JEDI_consensus.tsv'), 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow([
        'ASV', 'Domain_JEDI', 'Phylum_JEDI', 'Class_JEDI', 'Order_JEDI',
        'Family_JEDI', 'Genus_JEDI', 'Species_JEDI', 'Source'
    ])
    for asv in counts:
        silva_row = silva_tax.get(asv, {rank: '' for rank in RANKS})
        pr2_row = pr2_tax.get(asv, {rank: '' for rank in RANKS})
        domain, base, source = choose_assignment(silva_row, pr2_row)
        row = [asv, domain] + [clean_tax(base.get(rank, '')) for rank in RANKS[1:]] + [source]
        consensus.append(row)
        writer.writerow(row)
        total = sum(counts[asv].values())
        domain_summary[domain] += total
        for sample, value in counts[asv].items():
            domain_per_sample[sample][domain] += value

with open(os.path.join(integrated_dir, 'ASV_table_JEDI_counts.tsv'), 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['ASV'] + samples)
    for asv in counts:
        writer.writerow([asv] + [counts[asv][sample] for sample in samples])

with open(os.path.join(integrated_dir, 'ASV_table_JEDI_with_taxonomy.tsv'), 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow([
        'ASV', 'Domain_JEDI', 'Phylum_JEDI', 'Class_JEDI', 'Order_JEDI',
        'Family_JEDI', 'Genus_JEDI', 'Species_JEDI', 'Source'
    ] + samples + ['Sequence'])
    consensus_by_asv = {row[0]: row[1:] for row in consensus}
    for asv in counts:
        writer.writerow(
            [asv] + consensus_by_asv[asv] +
            [counts[asv][sample] for sample in samples] +
            [sequences.get(asv, '')]
        )

with open(os.path.join(integrated_dir, 'domain_summary.tsv'), 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['Domain', 'Total_reads'])
    for domain, total in sorted(domain_summary.items()):
        writer.writerow([domain, total])

all_domains = sorted({domain for values in domain_per_sample.values() for domain in values})
with open(os.path.join(integrated_dir, 'domain_counts_per_sample.tsv'), 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['Sample'] + all_domains)
    for sample in samples:
        writer.writerow([sample] + [domain_per_sample[sample].get(domain, 0)
                                    for domain in all_domains])

with open(os.path.join(integrated_dir, 'domain_relative_abundance_per_sample.tsv'), 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['Sample'] + all_domains)
    for sample in samples:
        total = sum(domain_per_sample[sample].values())
        values = [domain_per_sample[sample].get(domain, 0) / total if total else 0
                  for domain in all_domains]
        writer.writerow([sample] + values)

for domain in all_domains:
    safe_domain = re.sub(r'[^A-Za-z0-9_.-]+', '_', domain)
    path = os.path.join(integrated_dir, 'tables_by_domain', '{}_ASV_table.tsv'.format(safe_domain))
    domain_asvs = [row[0] for row in consensus if row[1] == domain]
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['ASV'] + samples)
        for asv in domain_asvs:
            writer.writerow([asv] + [counts[asv][sample] for sample in samples])

with open(os.path.join(integrated_dir, 'alpha_diversity_unrarefied.tsv'), 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['Sample', 'Observed_ASVs', 'Shannon', 'Simpson', 'Reads'])
    for sample in samples:
        values = [counts[asv][sample] for asv in counts if counts[asv][sample] > 0]
        reads = sum(values)
        observed = len(values)
        if reads:
            proportions = [value / reads for value in values]
            shannon = -sum(p * math.log(p) for p in proportions if p > 0)
            simpson = 1 - sum(p * p for p in proportions)
        else:
            shannon = 0.0
            simpson = 0.0
        writer.writerow([sample, observed, shannon, simpson, reads])

sample_totals = {sample: sum(counts[asv][sample] for asv in counts) for sample in samples}
positive_totals = [value for value in sample_totals.values() if value > 0]
max_depth = max(positive_totals) if positive_totals else 0
if max_depth:
    step = max(1000, max_depth // 20)
    depths = list(range(step, max_depth + 1, step))
    if max_depth not in depths:
        depths.append(max_depth)
else:
    depths = [0]

with open(os.path.join(integrated_dir, 'rarefaction_expected_ASVs.tsv'), 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['Sample', 'Depth', 'Expected_ASVs'])
    for sample in samples:
        total = sample_totals[sample]
        abundances = [counts[asv][sample] for asv in counts if counts[asv][sample] > 0]
        for depth in depths:
            if depth == 0 or total == 0:
                expected = 0.0
            elif depth >= total:
                expected = float(len(abundances))
            else:
                expected = 0.0
                denominator = (math.lgamma(total + 1) - math.lgamma(depth + 1) -
                               math.lgamma(total - depth + 1))
                for abundance in abundances:
                    if total - abundance >= depth:
                        numerator = (math.lgamma(total - abundance + 1) -
                                     math.lgamma(depth + 1) -
                                     math.lgamma(total - abundance - depth + 1))
                        p_zero = math.exp(numerator - denominator)
                    else:
                        p_zero = 0.0
                    expected += 1.0 - p_zero
            writer.writerow([sample, depth, expected])

print('Intégration terminée : {} ASV, {} échantillons'.format(len(counts), len(samples)))
PY_INTEGRATE
}

write_manifest() {
  log "Écriture du manifest final"
  cat > "${RUN_MANIFEST}" <<EOF_MANIFEST
script_version\t${SCRIPT_VERSION}
date\t$(date '+%F %T')
project_root\t${PROJECT_ROOT}
jedi_dir\t${JEDI_DIR}
raw_dir\t${RAW_DIR}
profile\t${PROFILE}
nfcore_version\t${NFCORE_VERSION}
forward_primer\t${FORWARD_PRIMER}
reverse_primer\t${REVERSE_PRIMER}
trunclenf\t${TRUNC_LEN_F}
trunclenr\t${TRUNC_LEN_R}
max_ee\t${MAX_EE}
trunc_qmin\t${TRUNC_QMIN}
trunc_rmin\t${TRUNC_RMIN}
silva_reference\t${SILVA_REF}
pr2_reference\t${PR2_REF_LABEL}
samplesheet\t${SAMPLESHEET}
silva_asv_fasta\t${SILVA_ASV_FASTA}
silva_asv_table\t${SILVA_ASV_TABLE}
silva_taxonomy\t${SILVA_TAXONOMY}
pr2_taxonomy\t${PR2_TAXONOMY_TSV}
integrated_dir\t${INTEGRATED_DIR}
EOF_MANIFEST
}

validate_final_outputs() {
  local required=(
    "${SILVA_ASV_FASTA}"
    "${SILVA_ASV_TABLE}"
    "${SILVA_TAXONOMY}"
    "${PR2_TAXONOMY_TSV}"
    "${INTEGRATED_DIR}/ASV_table_JEDI_counts.tsv"
    "${INTEGRATED_DIR}/ASV_table_JEDI_with_taxonomy.tsv"
    "${INTEGRATED_DIR}/taxonomy_JEDI_consensus.tsv"
    "${INTEGRATED_DIR}/taxonomy_SILVA_original.tsv"
    "${INTEGRATED_DIR}/taxonomy_PR2_original.tsv"
    "${INTEGRATED_DIR}/domain_summary.tsv"
    "${INTEGRATED_DIR}/domain_counts_per_sample.tsv"
    "${INTEGRATED_DIR}/domain_relative_abundance_per_sample.tsv"
    "${INTEGRATED_DIR}/alpha_diversity_unrarefied.tsv"
    "${INTEGRATED_DIR}/rarefaction_expected_ASVs.tsv"
    "${RUN_MANIFEST}"
  )

  local output
  for output in "${required[@]}"; do
    [[ -s "${output}" ]] || die "Sortie finale absente ou vide : ${output}"
  done

  printf 'SUCCESS\t%s\t%s\n' "$(date '+%F %T')" "${SCRIPT_VERSION}" > "${SUCCESS_FLAG}"
  if [[ "${KEEP_FAILURE_MARKER}" != "1" ]]; then
    rm -f "${FAIL_FLAG}" || true
  fi
  log "Validation finale OK"
}

main() {
  acquire_lock
  setup_environment
  activate_nextflow_java
  check_prerequisites
  write_nextflow_config
  preload_problematic_container
  build_samplesheet
  check_fastq_integrity

  log "Pipeline ${SCRIPT_VERSION}"
  log "Racine ${JEDI_DIR}"
  log "nf-core/ampliseq ${NFCORE_VERSION}"
  log "Profil ${PROFILE}"
  log "Classifieurs ${SILVA_REF} et ${PR2_REF_LABEL}"
  log "Amorces ${FORWARD_PRIMER} ${REVERSE_PRIMER}"

  run_silva_pipeline
  classify_pr2_from_silva_asvs
  integrate_outputs
  write_manifest
  validate_final_outputs

  log "Pipeline JEDI terminé avec succès"
  log "ASV SILVA : ${SILVA_ASV_FASTA}"
  log "Taxonomie SILVA : ${SILVA_TAXONOMY}"
  log "Taxonomie PR2 : ${PR2_TAXONOMY_TSV}"
  log "Table intégrée : ${INTEGRATED_DIR}/ASV_table_JEDI_with_taxonomy.tsv"
}

main "$@"
