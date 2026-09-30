#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

############################################
# JEDI cross-domain pipeline for Valormicro
# Complete fixed version
############################################

SCRIPT_VERSION="JEDI_VALORMICRO_FIXED_V3_2026-09-30"

PROJECT_ROOT="/nvme/bio/data_fungi/valormicro_V2"
RAW_DIR="${PROJECT_ROOT}/01_raw_data"
OLD_PIPE_DIR="${PROJECT_ROOT}/02_amplicon_pipeline"
OLD_SAMPLES_TSV="${OLD_PIPE_DIR}/04_database_files/samples.tsv"

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
LOG_FILE="${LOG_DIR}/jedi_pipeline_fixed_v3.log"

SAMPLESHEET="${INPUT_DIR}/samplesheet_jedi.tsv"
SAMPLESHEET_DEBUG="${INPUT_DIR}/samplesheet_jedi_debug.tsv"
SILVA_PARAMS="${INPUT_DIR}/params_silva_fixed_v3.yaml"
PR2_PARAMS="${INPUT_DIR}/params_pr2_fixed_v3.yaml"
NF_INFRA_CFG="${INPUT_DIR}/nextflow_infrastructure.config"

NFCORE_VERSION="2.18.0"
PROFILE="${PROFILE:-singularity}"
AMPLISEQ_REPO="${AMPLISEQ_REPO:-nf-core/ampliseq}"

PYTHON_BIN="${PYTHON_BIN:-python3}"

FORWARD_PRIMER="GTGYCAGCMGCCGCGGTAA"
REVERSE_PRIMER="CCGYCAATTYMTTTRAGTTT"

TRUNC_LEN_F="${TRUNC_LEN_F:-231}"
TRUNC_LEN_R="${TRUNC_LEN_R:-230}"

SILVA_REF="silva=138.2"
PR2_REF="pr2=5.0.0"

PRELOAD_IMAGE_DOCKER="docker://biocontainers/biocontainers:v1.2.0_cv1"
PRELOAD_IMAGE_SIF="${SINGULARITY_IMAGE_CACHE}/biocontainers_v1.2.0_cv1.sif"
PRELOAD_IMAGE_ALIAS1="${SINGULARITY_IMAGE_CACHE}/containers.biocontainers.pro-s3-SingImgsRepo-biocontainers-v1.2.0_cv1-biocontainers_v1.2.0_cv1.img.img"
PRELOAD_IMAGE_ALIAS2="${CONTAINER_ROOT}/containers.biocontainers.pro-s3-SingImgsRepo-biocontainers-v1.2.0_cv1-biocontainers_v1.2.0_cv1.img.img"

PULL_TIMEOUT="${PULL_TIMEOUT:-12 h}"

FORCE_SILVA="${FORCE_SILVA:-0}"
FORCE_PR2="${FORCE_PR2:-0}"
FORCE_ALL="${FORCE_ALL:-0}"
KEEP_FAILURE_MARKER="${KEEP_FAILURE_MARKER:-0}"

mkdir -p \
  "${INPUT_DIR}" \
  "${SILVA_DIR}" \
  "${PR2_DIR}" \
  "${INTEGRATED_DIR}" \
  "${REF_DIR}" \
  "${LOG_DIR}" \
  "${TMP_DIR}" \
  "${CONTAINER_ROOT}" \
  "${SINGULARITY_TMPDIR_LOCAL}" \
  "${SINGULARITY_LAYER_CACHE}" \
  "${SINGULARITY_IMAGE_CACHE}" \
  "${INTEGRATED_DIR}/tables_by_domain"

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
    else
        die "Un autre run JEDI semble déjà actif : ${LOCK_DIR}"
    fi
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
    export NXF_HOME="${JEDI_DIR}/.nextflow"
    export NXF_OFFLINE='false'
    mkdir -p "${NXF_HOME}" "${SINGULARITY_TMPDIR_LOCAL}" "${SINGULARITY_LAYER_CACHE}" "${SINGULARITY_IMAGE_CACHE}"
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
    require_cmd "${PYTHON_BIN}"
    require_cmd java

    [[ -d "${RAW_DIR}" ]] || die "Répertoire FASTQ absent : ${RAW_DIR}"
    [[ -f "${PROJECT_ROOT}/01_raw_data/00_infos_data.xlsx" || -f "${OLD_SAMPLES_TSV}" ]] || die \
      "Aucune source de métadonnées trouvée (ni ${OLD_SAMPLES_TSV}, ni 00_infos_data.xlsx)"
}

write_nextflow_config() {
    cat > "${NF_INFRA_CFG}" <<EOF
process {
  executor = 'local'
  scratch = false
  errorStrategy = 'terminate'
  maxRetries = 1
  cpus = 8
  memory = '64 GB'
  time = '72 h'
}

singularity {
  enabled = true
  autoMounts = true
  cacheDir = '${SINGULARITY_IMAGE_CACHE}'
  pullTimeout = '${PULL_TIMEOUT}'
}

apptainer {
  enabled = false
}

docker.enabled = false
podman.enabled = false
charliecloud.enabled = false

params {
  validate_params = true
}
EOF
}

singularity {
  enabled = true
  autoMounts = true
  cacheDir = '${SINGULARITY_IMAGE_CACHE}'
  pullTimeout = '${PULL_TIMEOUT}'
}

apptainer {
  enabled = true
  autoMounts = true
  cacheDir = '${SINGULARITY_IMAGE_CACHE}'
}

docker.enabled = false
podman.enabled = false
charliecloud.enabled = false

params {
  validate_params = true
}
EOF
}

preload_problematic_container() {
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

    "${PYTHON_BIN}" - "${OLD_SAMPLES_TSV}" "${RAW_DIR}" "${SAMPLESHEET}" "${SAMPLESHEET_DEBUG}" "${PROJECT_ROOT}/01_raw_data/00_infos_data.xlsx" <<'PY'
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
        value = "Sample"
    if not re.match(r'^[A-Za-z]', value):
        value = "S_" + value
    if not re.match(r'^[A-Za-z][A-Za-z0-9_]+$', value):
        raise SystemExit("ID nf-core invalide après normalisation : {}".format(value))
    return value

def resolve_fastq(raw_dir, filename):
    filename = os.path.basename(str(filename).strip())
    if not filename:
        return None
    direct = os.path.join(raw_dir, filename)
    if os.path.isfile(direct):
        return os.path.abspath(direct)

    found = []
    for root, _, files in os.walk(raw_dir):
        if filename in files:
            found.append(os.path.abspath(os.path.join(root, filename)))

    if len(found) == 1:
        return found[0]
    if len(found) == 0:
        return None
    raise SystemExit("FASTQ ambigu pour {}:\n{}".format(filename, "\n".join(found)))

def write_outputs(records, debug_records, out_tsv, debug_tsv):
    with open(out_tsv, "w", newline="") as oh:
        w = csv.writer(oh, delimiter="\t")
        w.writerow(["sample", "fastq_1", "fastq_2"])
        for rec in records:
            w.writerow([rec["sample"], rec["fastq_1"], rec["fastq_2"]])

    with open(debug_tsv, "w", newline="") as oh:
        w = csv.writer(oh, delimiter="\t")
        w.writerow(["sample", "fastq_1", "fastq_2", "source"])
        for rec in debug_records:
            w.writerow([rec["sample"], rec["fastq_1"], rec["fastq_2"], rec["source"]])

def build_from_old_tsv(path, raw_dir):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None

    records = []
    debug_records = []
    used = set()

    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit("En-tête absent dans {}".format(path))

        fields = {name.strip().lower(): name for name in reader.fieldnames if name}

        def field(*candidates):
            for candidate in candidates:
                if candidate.lower() in fields:
                    return fields[candidate.lower()]
            return None

        r1_col = field("R1")
        r2_col = field("R2")
        label_col = field("New_label", "Newlabel", "sample-id", "Label")

        missing = []
        if not r1_col:
            missing.append("R1")
        if not r2_col:
            missing.append("R2")
        if not label_col:
            missing.append("New_label/sample-id/Label")

        if missing:
            raise SystemExit("Colonnes requises absentes dans {} : {}".format(path, ", ".join(missing)))

        for line_num, row in enumerate(reader, start=2):
            r1_name = (row.get(r1_col) or "").strip()
            r2_name = (row.get(r2_col) or "").strip()
            label = (row.get(label_col) or "").strip()

            if not r1_name and not r2_name and not label:
                continue
            if not all([r1_name, r2_name, label]):
                raise SystemExit(
                    "Ligne incomplète dans {} à la ligne {} : R1='{}' R2='{}' label='{}'".format(
                        path, line_num, r1_name, r2_name, label
                    )
                )

            sample = nfcore_id(label)
            if sample in used:
                raise SystemExit("ID dupliqué après normalisation : {}".format(sample))

            r1_abs = resolve_fastq(raw_dir, r1_name)
            r2_abs = resolve_fastq(raw_dir, r2_name)
            if not r1_abs:
                raise SystemExit("FASTQ introuvable pour R1 '{}' (ligne {})".format(r1_name, line_num))
            if not r2_abs:
                raise SystemExit("FASTQ introuvable pour R2 '{}' (ligne {})".format(r2_name, line_num))

            used.add(sample)
            rec = {"sample": sample, "fastq_1": r1_abs, "fastq_2": r2_abs}
            records.append(rec)
            debug_records.append({**rec, "source": "samples.tsv"})

    return records, debug_records

def read_xlsx(path):
    ns = {'a': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}

    def col_to_idx(col):
        n = 0
        for c in col:
            if c.isalpha():
                n = n * 26 + (ord(c.upper()) - 64)
        return n - 1

    strings = []
    with zipfile.ZipFile(path) as z:
        if 'xl/sharedStrings.xml' in z.namelist():
            root = ET.fromstring(z.read('xl/sharedStrings.xml'))
            for si in root.findall('a:si', ns):
                texts = []
                for t in si.iterfind('.//a:t', ns):
                    texts.append(t.text or '')
                strings.append(''.join(texts))

        wb = ET.fromstring(z.read('xl/workbook.xml'))
        rel_root = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        rels = {r.attrib['Id']: r.attrib['Target'] for r in rel_root}

        sheets = []
        for s in wb.find('a:sheets', ns):
            rid = s.attrib.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
            target = rels[rid]
            if not target.startswith('xl/'):
                target = 'xl/' + target
            sheets.append((s.attrib['name'], target))

        all_rows = []
        for name, target in sheets:
            root = ET.fromstring(z.read(target))
            data = []
            sheetData = root.find('a:sheetData', ns)
            if sheetData is None:
                continue
            for row in sheetData.findall('a:row', ns):
                vals = {}
                for c in row.findall('a:c', ns):
                    ref = c.attrib.get('r', 'A1')
                    col = ''.join([x for x in ref if x.isalpha()])
                    idx = col_to_idx(col)
                    t = c.attrib.get('t')
                    v = c.find('a:v', ns)
                    val = ''
                    if v is not None and v.text is not None:
                        val = v.text
                        if t == 's':
                            val = strings[int(val)]
                    else:
                        isel = c.find('a:is', ns)
                        if isel is not None:
                            txt = [x.text or '' for x in isel.iterfind('.//a:t', ns)]
                            val = ''.join(txt)
                    vals[idx] = str(val).strip()

                if vals:
                    maxidx = max(vals)
                    rowvals = [''] * (maxidx + 1)
                    for i, v in vals.items():
                        rowvals[i] = v
                    data.append(rowvals)
            all_rows.append((name, data))
    return all_rows

def build_from_xlsx(xlsx, raw_dir):
    if not os.path.isfile(xlsx):
        raise SystemExit("Fichier Excel introuvable : {}".format(xlsx))

    def clean_header(x):
        return re.sub(r'\s+', ' ', str(x).strip().lower())

    def header_key(value):
        value = str(value).replace("\xa0", " ")
        value = value.strip().lower()
        return re.sub(r"[^a-z0-9]", "", value)

    rows = read_xlsx(xlsx)
    pairs = []
    debug_records = []
    seen = set()

    for sheet_name, sheet in rows:
        if not sheet:
            continue

        header = None
        header_idx = None

        for i, row in enumerate(sheet[:20]):
            norm_keys = [header_key(c) for c in row]
            if "r1" in norm_keys and "r2" in norm_keys:
                header = row
                header_idx = i
                break

        if header is None:
            continue

        header_keys = [header_key(x) for x in header]

        def find_col(accepted_names):
            accepted = {header_key(x) for x in accepted_names}
            for j, key in enumerate(header_keys):
                if key in accepted:
                    return j
            return None

        r1_col = find_col(("R1", "Read1", "Forward", "Forward read"))
        r2_col = find_col(("R2", "Read2", "Reverse", "Reverse read"))
        label_col = find_col((
            "Newlabel", "New label", "New_label",
            "Sample", "SampleID", "Sample ID",
            "Echantillon", "Échantillon", "Label"
        ))

        if r1_col is None or r2_col is None or label_col is None:
            continue

        for row_number, row in enumerate(sheet[header_idx + 1:], start=header_idx + 2):
            if not any(str(x).strip() for x in row):
                continue

            r1_name = row[r1_col].strip() if r1_col < len(row) else ""
            r2_name = row[r2_col].strip() if r2_col < len(row) else ""
            label = row[label_col].strip() if label_col < len(row) else ""

            if not r1_name and not r2_name and not label:
                continue

            if not all((r1_name, r2_name, label)):
                raise SystemExit(
                    "Feuille '{}', ligne Excel {} : R1, R2 et label doivent tous être renseignés.".format(
                        sheet_name, row_number
                    )
                )

            sample = nfcore_id(label)
            if sample in seen:
                raise SystemExit("ID nf-core dupliqué après normalisation : {}".format(sample))

            r1 = resolve_fastq(raw_dir, r1_name)
            r2 = resolve_fastq(raw_dir, r2_name)

            if not r1 or not r2:
                raise SystemExit(
                    "Feuille '{}', ligne Excel {} : FASTQ introuvable(s) pour R1='{}' R2='{}'".format(
                        sheet_name, row_number, r1_name, r2_name
                    )
                )

            seen.add(sample)
            rec = {"sample": sample, "fastq_1": r1, "fastq_2": r2}
            pairs.append(rec)
            debug_records.append({**rec, "source": "00_infos_data.xlsx::{}".format(sheet_name)})

    if not pairs:
        raise SystemExit("Aucune paire valide n'a été créée depuis {}".format(xlsx))

    return pairs, debug_records

records = None
debug_records = None

if os.path.isfile(old_tsv) and os.path.getsize(old_tsv) > 0:
    records, debug_records = build_from_old_tsv(old_tsv, raw_dir)
else:
    records, debug_records = build_from_xlsx(xlsx, raw_dir)

if not records:
    raise SystemExit("Aucun échantillon valide détecté.")

write_outputs(records, debug_records, out_tsv, debug_tsv)
print("OK : {} échantillons valides écrits.".format(len(records)))
PY

    sed -i 's/\r$//' "${SAMPLESHEET}" "${SAMPLESHEET_DEBUG}"

    local n
    n="$(awk 'END{print NR-1}' "${SAMPLESHEET}")"
    [[ "${n}" -gt 0 ]] || die "Samplesheet vide"
    log "${n} échantillons écrits dans ${SAMPLESHEET}"
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
    cat > "${SILVA_PARAMS}" <<EOF
input: "${SAMPLESHEET}"
outdir: "${SILVA_DIR}"
analysis: "dada2"
skip_cutadapt: false
cutadapt_minimum_length: 50
FW_primer: "${FORWARD_PRIMER}"
RV_primer: "${REVERSE_PRIMER}"
dada_ref_taxonomy: "silva"
dada_skip_exception: false
retain_untrimmed: false
trunc_qmin: 2
trunc_len_f: ${TRUNC_LEN_F}
trunc_len_r: ${TRUNC_LEN_R}
max_ee_f: 2
max_ee_r: 2
skip_dada_quality: false
dada2_min_fold_parent_over_abundance: 1.0
dada2_pooling_method: "independent"
merge_pairs: true
consensus_merge: true
pseudo_pool: false
sample_inference: false
qiime_metadata_category: ""
remove_taxids: ""
allow_multiple_taxonomies: true
skip_barplot: false
skip_abundance_tables: false
skip_qiime: false
skip_phylogeny: true
EOF
}

write_pr2_params() {
    cat > "${PR2_PARAMS}" <<EOF
input: "${SAMPLESHEET}"
outdir: "${PR2_DIR}"
analysis: "dada2"
skip_cutadapt: false
cutadapt_minimum_length: 50
FW_primer: "${FORWARD_PRIMER}"
RV_primer: "${REVERSE_PRIMER}"
dada_ref_taxonomy: "pr2"
dada_skip_exception: false
retain_untrimmed: false
trunc_qmin: 2
trunc_len_f: ${TRUNC_LEN_F}
trunc_len_r: ${TRUNC_LEN_R}
max_ee_f: 2
max_ee_r: 2
skip_dada_quality: false
dada2_min_fold_parent_over_abundance: 1.0
dada2_pooling_method: "independent"
merge_pairs: true
consensus_merge: true
pseudo_pool: false
sample_inference: false
qiime_metadata_category: ""
remove_taxids: ""
allow_multiple_taxonomies: true
skip_barplot: false
skip_abundance_tables: false
skip_qiime: false
skip_phylogeny: true
EOF
}

has_silva_success() {
    [[ -s "${SILVA_DIR}/dada2/table.tsv" ]] && \
    [[ -s "${SILVA_DIR}/dada2/representative_sequences.fasta" ]] && \
    find "${SILVA_DIR}" -type f \( -iname '*taxonomy*.tsv' -o -iname '*tax*.tsv' \) | grep -q .
}

has_pr2_success() {
    [[ -s "${PR2_DIR}/dada2/table.tsv" ]] && \
    [[ -s "${PR2_DIR}/dada2/representative_sequences.fasta" ]] && \
    find "${PR2_DIR}" -type f \( -iname '*taxonomy*.tsv' -o -iname '*tax*.tsv' \) | grep -q .
}

find_taxonomy_file() {
    local root="$1"
    find "${root}" -type f \( -iname '*taxonomy*.tsv' -o -iname '*tax*.tsv' \) \
        ! -iname '*barplot*' \
        ! -iname '*krona*' \
        | head -n 1
}

run_nfcore() {
    local params_file="$1"
    local work_subdir="$2"
    local label="$3"

    mkdir -p "$(dirname "${work_subdir}")" "${work_subdir}"

    log "Lancement nf-core/ampliseq ${label}"
    log "Paramètres : ${params_file}"

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

integrate_results() {
    log "Intégration SILVA + PR2"

    local silva_asv="${SILVA_DIR}/dada2/representative_sequences.fasta"
    local silva_table="${SILVA_DIR}/dada2/table.tsv"
    local silva_tax
    local pr2_tax

    silva_tax="$(find_taxonomy_file "${SILVA_DIR}")"
    pr2_tax="$(find_taxonomy_file "${PR2_DIR}")"

    [[ -s "${silva_asv}" ]] || die "FASTA ASV SILVA absent : ${silva_asv}"
    [[ -s "${silva_table}" ]] || die "Table ASV SILVA absente : ${silva_table}"
    [[ -n "${silva_tax}" && -s "${silva_tax}" ]] || die "Taxonomie SILVA introuvable"
    [[ -n "${pr2_tax}" && -s "${pr2_tax}" ]] || die "Taxonomie PR2 introuvable"

    "${PYTHON_BIN}" - "${silva_asv}" "${silva_tax}" "${pr2_tax}" "${silva_table}" "${INTEGRATED_DIR}" <<'PY'
import csv
import math
import os
import re
import sys
from collections import OrderedDict, defaultdict

silva_asv, silva_tax, pr2_tax, silva_table, integrated_dir = sys.argv[1:6]
os.makedirs(os.path.join(integrated_dir, "tables_by_domain"), exist_ok=True)

def parse_fasta(path):
    seqs = OrderedDict()
    cur = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = ''.join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if cur is not None:
        seqs[cur] = ''.join(buf)
    return seqs

def sniff_delim(header_line):
    return '\t' if header_line.count('\t') >= header_line.count(',') else ','

def normalize_rank_name(x):
    x = (x or '').strip().lower()
    x = x.replace('kingdom', 'domain')
    return x

def parse_tax_table(path):
    with open(path) as fh:
        first = fh.readline()
        if not first:
            return {}
    delim = sniff_delim(first)

    data = {}
    with open(path) as fh:
        reader = csv.reader(fh, delimiter=delim)
        rows = list(reader)

    if not rows:
        return data

    header = [normalize_rank_name(x) for x in rows[0]]
    idx_asv = None
    for i, h in enumerate(header):
        if h in ('featureid', 'feature id', 'asv', 'asv_id', 'sequence', 'otu', '#otuid', 'id'):
            idx_asv = i
            break
    if idx_asv is None:
        idx_asv = 0

    rank_names = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    rank_idx = {}
    for r in rank_names:
        for i, h in enumerate(header):
            if h == r or h.endswith('_' + r) or h.startswith(r + '_'):
                rank_idx[r] = i
                break

    tax_idx = None
    for i, h in enumerate(header):
        if 'taxonomy' in h or h in ('tax', 'taxon'):
            tax_idx = i
            break

    for row in rows[1:]:
        if not row or idx_asv >= len(row):
            continue
        asv = row[idx_asv].strip()
        if not asv:
            continue

        ranks = {r: '' for r in rank_names}
        if rank_idx:
            for r, i in rank_idx.items():
                if i < len(row):
                    ranks[r] = row[i].strip()
        elif tax_idx is not None and tax_idx < len(row):
            tax = row[tax_idx].strip()
            parts = re.split(r'[;,]\s*', tax)
            cleaned = []
            for p in parts:
                p = re.sub(r'^[dkpcofgs]__', '', p)
                cleaned.append(p)
            for i, r in enumerate(rank_names):
                if i < len(cleaned):
                    ranks[r] = cleaned[i]

        data[asv] = ranks

    return data

def parse_count_table(path):
    with open(path) as fh:
        first = fh.readline().rstrip('\n')
    delim = sniff_delim(first)

    counts = OrderedDict()
    with open(path) as fh:
        reader = csv.reader(fh, delimiter=delim)
        header = next(reader)
        samples = [x.strip() for x in header[1:]]
        for row in reader:
            if not row:
                continue
            asv = row[0].strip()
            vals = []
            for x in row[1:1 + len(samples)]:
                x = x.strip()
                vals.append(int(float(x)) if x not in ('', 'NA', 'nan') else 0)
            counts[asv] = dict(zip(samples, vals))
    return samples, counts

def clean_tax(x):
    x = (x or '').strip()
    if x.lower() in ('', 'na', 'nan', 'none', 'unassigned', 'unknown'):
        return ''
    return x

def choose_domain(silva, pr2):
    sdom = clean_tax(silva.get('domain', ''))
    pdom = clean_tax(pr2.get('domain', ''))
    text = ' '.join([clean_tax(v) for v in silva.values()]).lower()

    if 'mitochond' in text:
        return 'Mitochondria'
    if 'chloroplast' in text or 'plastid' in text:
        return 'Plastid'
    if sdom.lower() == 'bacteria':
        return 'Bacteria'
    if sdom.lower() == 'archaea':
        return 'Archaea'
    if pdom:
        return 'Eukaryota'
    if sdom and sdom.lower() not in ('eukaryota', 'eukaryote', 'opisthokonta'):
        return sdom
    return 'Unresolved'

seqs = parse_fasta(silva_asv)
silva = parse_tax_table(silva_tax)
pr2 = parse_tax_table(pr2_tax)
samples, counts = parse_count_table(silva_table)

consensus = []
domain_summary = defaultdict(int)
domain_per_sample = {s: defaultdict(int) for s in samples}

with open(os.path.join(integrated_dir, 'taxonomy_SILVA_original.tsv'), 'w', newline='') as oh:
    w = csv.writer(oh, delimiter='\t')
    w.writerow(['ASV', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
    for asv in counts:
        r = silva.get(asv, {})
        w.writerow([asv] + [clean_tax(r.get(k, '')) for k in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']])

with open(os.path.join(integrated_dir, 'taxonomy_PR2_original.tsv'), 'w', newline='') as oh:
    w = csv.writer(oh, delimiter='\t')
    w.writerow(['ASV', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
    for asv in counts:
        r = pr2.get(asv, {})
        w.writerow([asv] + [clean_tax(r.get(k, '')) for k in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']])

with open(os.path.join(integrated_dir, 'taxonomy_JEDI_consensus.tsv'), 'w', newline='') as oh:
    w = csv.writer(oh, delimiter='\t')
    w.writerow(['ASV', 'Domain_JEDI', 'Phylum_JEDI', 'Class_JEDI', 'Order_JEDI', 'Family_JEDI', 'Genus_JEDI', 'Species_JEDI', 'Source'])
    for asv in counts:
        s = silva.get(asv, {k: '' for k in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']})
        p = pr2.get(asv, {k: '' for k in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']})
        domain = choose_domain(s, p)

        if domain in ('Bacteria', 'Archaea', 'Mitochondria', 'Plastid'):
            src = 'SILVA'
            base = s
        elif domain == 'Eukaryota':
            src = 'PR2'
            base = p
        else:
            src = 'UNRESOLVED'
            base = s if any(clean_tax(v) for v in s.values()) else p

        row = [asv, domain] + [clean_tax(base.get(k, '')) for k in ['phylum', 'class', 'order', 'family', 'genus', 'species']] + [src]
        consensus.append(row)
        w.writerow(row)

        total = sum(counts[asv].values())
        domain_summary[domain] += total
        for sample, val in counts[asv].items():
            domain_per_sample[sample][domain] += val

with open(os.path.join(integrated_dir, 'ASV_table_JEDI_counts.tsv'), 'w', newline='') as oh:
    w = csv.writer(oh, delimiter='\t')
    w.writerow(['ASV'] + samples)
    for asv in counts:
        w.writerow([asv] + [counts[asv][s] for s in samples])

with open(os.path.join(integrated_dir, 'ASV_table_JEDI_with_taxonomy.tsv'), 'w', newline='') as oh:
    w = csv.writer(oh, delimiter='\t')
    w.writerow(['ASV', 'Domain_JEDI', 'Phylum_JEDI', 'Class_JEDI', 'Order_JEDI', 'Family_JEDI', 'Genus_JEDI', 'Species_JEDI', 'Source'] + samples + ['Sequence'])
    cdict = {r[0]: r[1:] for r in consensus}
    for asv in counts:
        row = cdict[asv]
        w.writerow([asv] + row + [counts[asv][s] for s in samples] + [seqs.get(asv, '')])

with open(os.path.join(integrated_dir, 'domain_summary.tsv'), 'w', newline='') as oh:
    w = csv.writer(oh, delimiter='\t')
    w.writerow(['Domain', 'Total_reads'])
    for dom, total in sorted(domain_summary.items()):
        w.writerow([dom, total])

all_domains = sorted({d for smp in domain_per_sample.values() for d in smp.keys()})

with open(os.path.join(integrated_dir, 'domain_counts_per_sample.tsv'), 'w', newline='') as oh:
    w = csv.writer(oh, delimiter='\t')
    w.writerow(['Sample'] + all_domains)
    for s in samples:
        w.writerow([s] + [domain_per_sample[s].get(d, 0) for d in all_domains])

with open(os.path.join(integrated_dir, 'domain_relative_abundance_per_sample.tsv'), 'w', newline='') as oh:
    w = csv.writer(oh, delimiter='\t')
    w.writerow(['Sample'] + all_domains)
    for s in samples:
        total = sum(domain_per_sample[s].values())
        vals = [(domain_per_sample[s].get(d, 0) / total if total else 0) for d in all_domains]
        w.writerow([s] + vals)

for dom in all_domains:
    path = os.path.join(integrated_dir, 'tables_by_domain', f'{dom}_ASV_table.tsv')
    dom_asvs = [r[0] for r in consensus if r[1] == dom]
    with open(path, 'w', newline='') as oh:
        w = csv.writer(oh, delimiter='\t')
        w.writerow(['ASV'] + samples)
        for asv in dom_asvs:
            w.writerow([asv] + [counts[asv][s] for s in samples])

with open(os.path.join(integrated_dir, 'alpha_diversity_unrarefied.tsv'), 'w', newline='') as oh:
    w = csv.writer(oh, delimiter='\t')
    w.writerow(['Sample', 'Observed_ASVs', 'Shannon', 'Simpson', 'Reads'])
    for s in samples:
        vec = [counts[a][s] for a in counts if counts[a][s] > 0]
        reads = sum(vec)
        observed = len(vec)
        if reads:
            ps = [v / reads for v in vec]
            sh = -sum(p * math.log(p) for p in ps if p > 0)
            sim = 1 - sum(p * p for p in ps)
        else:
            sh = 0.0
            sim = 0.0
        w.writerow([s, observed, sh, sim, reads])

with open(os.path.join(integrated_dir, 'rarefaction_expected_ASVs.tsv'), 'w', newline='') as oh:
    w = csv.writer(oh, delimiter='\t')
    totals = {s: sum(counts[a][s] for a in counts) for s in samples}
    max_depth = max(totals.values()) if totals else 0
    if max_depth == 0:
        depths = [0]
    else:
        step = max(1000, max_depth // 20)
        depths = list(range(step, max_depth + 1, step))
    w.writerow(['Sample', 'Depth', 'Expected_ASVs'])
    for s in samples:
        total = totals[s]
        abund = [counts[a][s] for a in counts if counts[a][s] > 0]
        for m in depths:
            if total == 0 or m == 0:
                exp = 0.0
            elif m >= total:
                exp = float(len(abund))
            else:
                exp = 0.0
                for n in abund:
                    if total - n >= m:
                        num = math.lgamma(total - n + 1) - math.lgamma(m + 1) - math.lgamma(total - n - m + 1)
                        den = math.lgamma(total + 1) - math.lgamma(m + 1) - math.lgamma(total - m + 1)
                        p0 = math.exp(num - den)
                    else:
                        p0 = 0.0
                    exp += 1 - p0
            w.writerow([s, m, exp])
PY
}

write_manifest() {
    log "Écriture du manifest final"
    cat > "${RUN_MANIFEST}" <<EOF
script_version	${SCRIPT_VERSION}
date	$(date '+%F %T')
project_root	${PROJECT_ROOT}
samplesheet	${SAMPLESHEET}
silva_dir	${SILVA_DIR}
pr2_dir	${PR2_DIR}
integrated_dir	${INTEGRATED_DIR}
nfcore_repo	${AMPLISEQ_REPO}
nfcore_version	${NFCORE_VERSION}
profile	${PROFILE}
forward_primer	${FORWARD_PRIMER}
reverse_primer	${REVERSE_PRIMER}
trunc_len_f	${TRUNC_LEN_F}
trunc_len_r	${TRUNC_LEN_R}
EOF
}

mark_success() {
    cat > "${SUCCESS_FLAG}" <<EOF
SUCCESS	$(date '+%F %T')
SCRIPT_VERSION	${SCRIPT_VERSION}
EOF

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
    log "Classifieurs ${SILVA_REF} et ${PR2_REF}"
    log "Amorces ${FORWARD_PRIMER} / ${REVERSE_PRIMER}"

    if [[ "${FORCE_ALL}" == "1" ]]; then
        FORCE_SILVA=1
        FORCE_PR2=1
    fi

    if [[ "${FORCE_SILVA}" == "1" || ! $(has_silva_success; echo $?) -eq 0 ]]; then
        log "ÉTAPE 1/3 : DADA2/JEDI et taxonomie SILVA"
        write_silva_params
        run_nfcore "${SILVA_PARAMS}" "${JEDI_DIR}/work/silva" "SILVA"
    else
        log "ÉTAPE 1/3 : SILVA déjà disponible, réutilisation"
    fi

    has_silva_success || die "Les sorties SILVA minimales ne sont pas présentes après exécution"

    if [[ "${FORCE_PR2}" == "1" || ! $(has_pr2_success; echo $?) -eq 0 ]]; then
        log "ÉTAPE 2/3 : taxonomie PR2"
        write_pr2_params
        run_nfcore "${PR2_PARAMS}" "${JEDI_DIR}/work/pr2" "PR2"
    else
        log "ÉTAPE 2/3 : PR2 déjà disponible, réutilisation"
    fi

    has_pr2_success || die "Les sorties PR2 minimales ne sont pas présentes après exécution"

    log "ÉTAPE 3/3 : intégration JEDI"
    integrate_results
    write_manifest
    mark_success
    log "Pipeline JEDI terminé avec succès"
}

main "$@"
