#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

############################################
# JEDI cross-domain pipeline for Valormicro
# Fixed v2: robust Singularity cache + preloaded problematic image
############################################

SCRIPT_VERSION="JEDI_VALORMICRO_FIXED_V2_2026-09-30"
PROJECT_ROOT="/nvme/bio/data_fungi/valormicro_V2"
RAW_DIR="${PROJECT_ROOT}/01_raw_data"
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
LOG_FILE="${LOG_DIR}/jedi_pipeline_fixed_v2.log"
SAMPLESHEET="${INPUT_DIR}/samplesheet_jedi.tsv"
SILVA_PARAMS="${INPUT_DIR}/params_silva_fixed_v2.yaml"
PR2_PARAMS="${INPUT_DIR}/params_pr2_fixed_v2.yaml"
NF_INFRA_CFG="${INPUT_DIR}/nextflow_infrastructure.config"
NFCORE_VERSION="2.18.0"
PROFILE="singularity"
AMPLISEQ_REPO="nf-core/ampliseq"
PYTHON_BIN="python3"
FORWARD_PRIMER="GTGYCAGCMGCCGCGGTAA"
REVERSE_PRIMER="CCGYCAATTYMTTTRAGTTT"
TRUNC_LEN_F="231"
TRUNC_LEN_R="230"
SILVA_REF="silva=138.2"
PR2_REF="pr2=5.0.0"
PRELOAD_IMAGE_DOCKER="docker://biocontainers/biocontainers:v1.2.0_cv1"
PRELOAD_IMAGE_SIF="${SINGULARITY_IMAGE_CACHE}/biocontainers_v1.2.0_cv1.sif"
PRELOAD_IMAGE_ALIAS1="${SINGULARITY_IMAGE_CACHE}/containers.biocontainers.pro-s3-SingImgsRepo-biocontainers-v1.2.0_cv1-biocontainers_v1.2.0_cv1.img.img"
PRELOAD_IMAGE_ALIAS2="${CONTAINER_ROOT}/containers.biocontainers.pro-s3-SingImgsRepo-biocontainers-v1.2.0_cv1-biocontainers_v1.2.0_cv1.img.img"
PULL_TIMEOUT="12 h"
FORCE_SILVA="${FORCE_SILVA:-0}"
FORCE_PR2="${FORCE_PR2:-0}"
FORCE_ALL="${FORCE_ALL:-0}"
KEEP_FAILURE_MARKER="${KEEP_FAILURE_MARKER:-0}"

mkdir -p "${INPUT_DIR}" "${SILVA_DIR}" "${PR2_DIR}" "${INTEGRATED_DIR}" "${REF_DIR}" "${LOG_DIR}" "${TMP_DIR}" \
         "${CONTAINER_ROOT}" "${SINGULARITY_TMPDIR_LOCAL}" "${SINGULARITY_LAYER_CACHE}" "${SINGULARITY_IMAGE_CACHE}"

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
  printf 'FAILED\t%s\tline=%s\tcmd=%s\texit=%s\n' "$(date '+%F %T')" "${line_no}" "${cmd}" "${exit_code}" > "${FAIL_FLAG}" || true
  exit "${exit_code}"
}
trap 'on_error ${LINENO} "${BASH_COMMAND}"' ERR

cleanup_lock() {
  rm -rf "${LOCK_DIR}" || true
}
trap cleanup_lock EXIT

acquire_lock() {
  if mkdir "${LOCK_DIR}" 2>/dev/null; then
    printf '%s\n' "$${PPID}" > "${LOCK_DIR}/pid"
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
      conda activate nextflow-java
    fi
  fi

  export JAVA_HOME="$(dirname "$(dirname "$(readlink -f "$(command -v java)")")")"
  export NXF_JAVA_HOME="${JAVA_HOME}"
  unset JAVA_CMD
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
  [[ -f "${PROJECT_ROOT}/01_raw_data/00_infos_data.xlsx" ]] || die "Fichier metadata absent : ${PROJECT_ROOT}/01_raw_data/00_infos_data.xlsx"
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

apptainer {
  enabled = false
}

docker {
  enabled = false
}

cleanup = false
report {
  enabled = false
}
timeline {
  enabled = false
}
trace {
  enabled = true
}
dag {
  enabled = false
}
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
  log "Construction et validation du samplesheet depuis ${PROJECT_ROOT}/01_raw_data/00_infos_data.xlsx"
  "${PYTHON_BIN}" <<'PY'
import os, re, sys, json, csv, zipfile, xml.etree.ElementTree as ET
project_root = "/nvme/bio/data_fungi/valormicro_V2"
raw_dir = os.path.join(project_root, "01_raw_data")
outfile = os.path.join(project_root, "03_JEDI_pipeline", "00_inputs", "samplesheet_jedi.tsv")
xlsx = os.path.join(raw_dir, "00_infos_data.xlsx")
shared = os.path.join(project_root, "03_JEDI_pipeline", "00_inputs", "samplesheet_jedi_debug.tsv")
ns = {'a':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}

def col_to_idx(col):
    n=0
    for c in col:
        if c.isalpha():
            n=n*26+(ord(c.upper())-64)
    return n-1

def read_xlsx(path):
    with zipfile.ZipFile(path) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            root=ET.fromstring(z.read('xl/sharedStrings.xml'))
            for si in root.findall('a:si', ns):
                texts=[]
                for t in si.iterfind('.//a:t', ns):
                    texts.append(t.text or '')
                strings.append(''.join(texts))
        wb=ET.fromstring(z.read('xl/workbook.xml'))
        rel_root=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        rels={r.attrib['Id']:r.attrib['Target'] for r in rel_root}
        sheets=[]
        for s in wb.find('a:sheets', ns):
            rid=s.attrib.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
            target=rels[rid]
            if not target.startswith('xl/'):
                target='xl/'+target
            sheets.append((s.attrib['name'], target))
        all_rows=[]
        for name,target in sheets:
            root=ET.fromstring(z.read(target))
            data=[]
            sheetData=root.find('a:sheetData', ns)
            if sheetData is None:
                continue
            for row in sheetData.findall('a:row', ns):
                vals={}
                for c in row.findall('a:c', ns):
                    ref=c.attrib.get('r','A1')
                    col=''.join([x for x in ref if x.isalpha()])
                    idx=col_to_idx(col)
                    t=c.attrib.get('t')
                    v=c.find('a:v', ns)
                    val=''
                    if v is not None and v.text is not None:
                        val=v.text
                        if t=='s':
                            val=strings[int(val)]
                    else:
                        isel=c.find('a:is', ns)
                        if isel is not None:
                            txt=[x.text or '' for x in isel.iterfind('.//a:t', ns)]
                            val=''.join(txt)
                    vals[idx]=val
                if vals:
                    maxidx=max(vals)
                    rowvals=['']*(maxidx+1)
                    for i,v in vals.items():
                        rowvals[i]=str(v).strip()
                    data.append(rowvals)
            all_rows.append((name,data))
        return all_rows

fastq_r1 = {}
fastq_r2 = {}
for fn in os.listdir(raw_dir):
    if fn.endswith('.fastq.gz') or fn.endswith('.fq.gz'):
        full=os.path.join(raw_dir, fn)
        m1=re.match(r'(.+?)(_R?1(?:_001)?)(\.f(?:ast)?q\.gz)$', fn)
        m2=re.match(r'(.+?)(_R?2(?:_001)?)(\.f(?:ast)?q\.gz)$', fn)
        if m1:
            fastq_r1[m1.group(1)] = full
        elif m2:
            fastq_r2[m2.group(1)] = full

pairs=[]
for key in sorted(set(fastq_r1) & set(fastq_r2)):
    sample=key
    sample=re.sub(r'[^A-Za-z0-9_.-]+', '_', sample)
    pairs.append((sample, fastq_r1[key], fastq_r2[key]))

if not pairs:
    raise SystemExit('Aucune paire R1/R2 détectée dans 01_raw_data')

rows = read_xlsx(xlsx)
negative_tokens = {'blank','neg','negative','ntc','control','controle','ctrl'}
metadata = {}
for _, sheet in rows:
    if not sheet:
        continue
    header=None
    for row in sheet[:10]:
        norm=[re.sub(r'\s+',' ',c.strip().lower()) for c in row]
        if any('sample' in x or 'echant' in x or 'échant' in x for x in norm):
            header=norm
            break
    if header is None:
        continue
    header_idx=sheet.index(row)
    for r in sheet[header_idx+1:]:
        if not any(str(x).strip() for x in r):
            continue
        rec={header[i]: r[i].strip() if i < len(r) else '' for i in range(len(header))}
        vals=' '.join(rec.values()).lower()
        sid=None
        for k,v in rec.items():
            if 'sample' in k or 'echant' in k or 'échant' in k:
                sid=v.strip()
                break
        if sid:
            sid_clean=re.sub(r'[^A-Za-z0-9_.-]+', '_', sid)
            metadata[sid_clean]={'negative':'yes' if any(tok in vals for tok in negative_tokens) else 'no'}

with open(outfile, 'w', newline='') as oh, open(shared,'w',newline='') as dbg:
    w=csv.writer(oh, delimiter='\t')
    d=csv.writer(dbg, delimiter='\t')
    w.writerow(['sample','fastq_1','fastq_2'])
    d.writerow(['sample','fastq_1','fastq_2','negative'])
    for sample,r1,r2 in pairs:
        w.writerow([sample,r1,r2])
        d.writerow([sample,r1,r2,metadata.get(sample,{}).get('negative','no')])
print(len(pairs))
PY
  local n
  n=$(awk 'END{print NR-1}' "${SAMPLESHEET}")
  [[ "${n}" -gt 0 ]] || die "Samplesheet vide"
  log "${n} échantillons écrits dans ${SAMPLESHEET}"
  sed -i 's/\r$//' "${SAMPLESHEET}"
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

    r1="${r1%$'\r'}"
    r2="${r2%$'\r'}"

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
trim_left_f: 0
trim_left_r: 0
trunc_len_f: ${TRUNC_LEN_F}
trunc_len_r: ${TRUNC_LEN_R}
max_ee_f: 2
max_ee_r: 2
trunc_q: 2
skip_fastqc: false
with_cutadapt: true
cutadapt_min_overlap: 5
cutadapt_error_rate: 0.1
denoise: "dada2"
single_end: false
multiple_sequencing_runs: false
pool_dada2: false
mergepairs_strategy: "consensus"
dada_ref_taxonomy: "${SILVA_REF}"
cut_dada_ref_taxonomy: true
exclude_taxa: "none"
skip_barrnap: false
skip_qiime: true
skip_qiime_downstream: true
skip_dada_addspecies: true
skip_abundance_tables: true
skip_alpha_rarefaction: true
skip_diversity_indices: true
ancombc_formula: "1"
report_title: "JEDI valormicro - ASV et SILVA"
save_intermediates: true
trunclenf: ${TRUNC_LEN_F}
trunclenr: ${TRUNC_LEN_R}
EOF_SILVA
}

write_pr2_params() {
  local asv_fasta="${SILVA_DIR}/dada2/ASV_seqs.fasta"
  [[ -s "${asv_fasta}" ]] || die "ASV fasta absent pour PR2 : ${asv_fasta}"
  cat > "${PR2_PARAMS}" <<EOF_PR2
input: "${asv_fasta}"
outdir: "${PR2_DIR}"
denoise: false
single_end: true
dada_ref_taxonomy: "${PR2_REF}"
cut_dada_ref_taxonomy: false
skip_barrnap: true
skip_fastqc: true
skip_qiime: true
skip_qiime_downstream: true
skip_dada_addspecies: true
skip_abundance_tables: true
skip_alpha_rarefaction: true
skip_diversity_indices: true
save_intermediates: true
report_title: "JEDI valormicro - PR2 sur ASV pré-calculés"
EOF_PR2
}

run_nfcore() {
  local params_file="$1"
  local work_subdir="$2"
  local label="$3"
  log "Lancement ${AMPLISEQ_REPO} ${label}"
  log "Paramètres : ${params_file}"
  log "Work : ${work_subdir}"
  ( cd "${JEDI_DIR}"; nextflow run "${AMPLISEQ_REPO}" -r "${NFCORE_VERSION}" -profile "${PROFILE}" -params-file "${params_file}" -work-dir "${work_subdir}" -c "${NF_INFRA_CFG}" -resume -ansi-log false )
}

has_silva_success() {
  local asv_fasta="${SILVA_DIR}/dada2/ASV_seqs.fasta"
  local asv_table="${SILVA_DIR}/dada2/table.tsv"
  local taxonomy_file=""

  [[ -s "${asv_fasta}" ]] || return 1
  [[ -s "${asv_table}" ]] || return 1

  taxonomy_file="$(
    find "${SILVA_DIR}" -type f \
      \( -iname '*taxonomy*.tsv' -o -iname '*tax*.tsv' -o -iname '*classification*.tsv' \) \
      ! -size 0c \
      | head -n 1
  )"

  [[ -n "${taxonomy_file}" ]]
}

has_pr2_success() {
  [[ -d "${PR2_DIR}" ]] && find "${PR2_DIR}" -type f | grep -q .
}

locate_taxonomy_file() {
  local root="$1"
  find "${root}" -type f \( -iname '*taxonomy*.tsv' -o -iname '*tax*.tsv' -o -iname '*classification*.tsv' \) | head -n 1
}

integrate_outputs() {
  log "Intégration SILVA + PR2"
  local silva_table="${SILVA_DIR}/dada2/table.tsv"
  local silva_asv="${SILVA_DIR}/dada2/ASV_seqs.fasta"
  [[ -s "${silva_table}" ]] || die "Table DADA2 absente : ${silva_table}"
  [[ -s "${silva_asv}" ]] || die "FASTA ASV absente : ${silva_asv}"
  local silva_tax
  silva_tax="$(locate_taxonomy_file "${SILVA_DIR}")"
  [[ -n "${silva_tax}" ]] || die "Fichier taxonomie SILVA introuvable"
  local pr2_tax
  pr2_tax="$(locate_taxonomy_file "${PR2_DIR}")"
  [[ -n "${pr2_tax}" ]] || die "Fichier taxonomie PR2 introuvable"

  mkdir -p "${INTEGRATED_DIR}/tables_by_domain"

  SILVA_TABLE="${silva_table}" SILVA_ASV="${silva_asv}" SILVA_TAX="${silva_tax}" PR2_TAX="${pr2_tax}" INTEGRATED_DIR="${INTEGRATED_DIR}" "${PYTHON_BIN}" <<'PY'
import os, re, math, csv
from collections import defaultdict, OrderedDict

integrated_dir = os.environ['INTEGRATED_DIR']
silva_table = os.environ['SILVA_TABLE']
silva_asv = os.environ['SILVA_ASV']
silva_tax = os.environ['SILVA_TAX']
pr2_tax = os.environ['PR2_TAX']

os.makedirs(integrated_dir, exist_ok=True)
os.makedirs(os.path.join(integrated_dir, 'tables_by_domain'), exist_ok=True)


def parse_fasta(path):
    seqs={}
    cur=None
    buf=[]
    with open(path) as fh:
        for line in fh:
            line=line.rstrip('\n')
            if not line: continue
            if line.startswith('>'):
                if cur is not None:
                    seqs[cur]=''.join(buf)
                cur=line[1:].split()[0]
                buf=[]
            else:
                buf.append(line)
        if cur is not None:
            seqs[cur]=''.join(buf)
    return seqs


def sniff_delim(header_line):
    return '\t' if header_line.count('\t') >= header_line.count(',') else ','


def normalize_rank_name(x):
    x=(x or '').strip().lower()
    x=x.replace('taxon','taxonomy').replace('kingdom','domain')
    return x


def parse_tax_table(path):
    with open(path) as fh:
        first=fh.readline()
        if not first:
            return {}
        delim=sniff_delim(first)
    data={}
    with open(path) as fh:
        reader=csv.reader(fh, delimiter=delim)
        rows=list(reader)
    if not rows:
        return data
    header=[normalize_rank_name(x) for x in rows[0]]
    idx_asv=None
    for i,h in enumerate(header):
        if h in ('featureid','feature id','asv','asv_id','sequence','otu','#otuid','id'):
            idx_asv=i
            break
    if idx_asv is None:
        idx_asv=0
    rank_names=['domain','phylum','class','order','family','genus','species']
    rank_idx={}
    for r in rank_names:
        for i,h in enumerate(header):
            if h==r or h.endswith('_'+r) or h.startswith(r+'_'):
                rank_idx[r]=i
                break
    tax_idx=None
    for i,h in enumerate(header):
        if 'taxonomy' in h or h in ('tax','taxon'):
            tax_idx=i
            break
    for row in rows[1:]:
        if not row or idx_asv >= len(row):
            continue
        asv=row[idx_asv].strip()
        if not asv:
            continue
        ranks={r:'' for r in rank_names}
        if rank_idx:
            for r,i in rank_idx.items():
                if i < len(row):
                    ranks[r]=row[i].strip()
        elif tax_idx is not None and tax_idx < len(row):
            tax=row[tax_idx].strip()
            parts=re.split(r'[;,]\s*', tax)
            cleaned=[]
            for p in parts:
                p=re.sub(r'^[dkpcofgs]__','',p)
                cleaned.append(p)
            for i,r in enumerate(rank_names):
                if i < len(cleaned):
                    ranks[r]=cleaned[i]
        data[asv]=ranks
    return data


def parse_count_table(path):
    with open(path) as fh:
        first=fh.readline().rstrip('\n')
        delim=sniff_delim(first)
    counts=OrderedDict()
    with open(path) as fh:
        reader=csv.reader(fh, delimiter=delim)
        header=next(reader)
        samples=[x.strip() for x in header[1:]]
        for row in reader:
            if not row:
                continue
            asv=row[0].strip()
            vals=[]
            for x in row[1:1+len(samples)]:
                x=x.strip()
                vals.append(int(float(x)) if x not in ('','NA','nan') else 0)
            counts[asv]=dict(zip(samples, vals))
    return samples, counts


def clean_tax(x):
    x=(x or '').strip()
    if x.lower() in ('', 'na', 'nan', 'none', 'unassigned', 'unknown'):
        return ''
    return x


def choose_domain(silva, pr2):
    sdom=clean_tax(silva.get('domain',''))
    pdom=clean_tax(pr2.get('domain',''))
    text=' '.join([clean_tax(v) for v in silva.values()]).lower()
    if 'mitochond' in text:
        return 'Mitochondria'
    if 'chloroplast' in text or 'plastid' in text:
        return 'Plastid'
    if sdom.lower() in ('bacteria','archaea'):
        return sdom.capitalize() if sdom.lower()=='bacteria' else 'Archaea'
    if pdom:
        return 'Eukaryota'
    if sdom and sdom.lower() not in ('eukaryota','eukaryote','opisthokonta'):
        return sdom
    return 'Unresolved'

seqs=parse_fasta(silva_asv)
silva=parse_tax_table(silva_tax)
pr2=parse_tax_table(pr2_tax)
samples, counts=parse_count_table(silva_table)

consensus=[]
domain_summary=defaultdict(int)
domain_per_sample={s:defaultdict(int) for s in samples}

with open(os.path.join(integrated_dir, 'taxonomy_SILVA_original.tsv'), 'w', newline='') as oh:
    w=csv.writer(oh, delimiter='\t')
    w.writerow(['ASV','domain','phylum','class','order','family','genus','species'])
    for asv in counts:
        r=silva.get(asv,{})
        w.writerow([asv]+[clean_tax(r.get(k,'')) for k in ['domain','phylum','class','order','family','genus','species']])

with open(os.path.join(integrated_dir, 'taxonomy_PR2_original.tsv'), 'w', newline='') as oh:
    w=csv.writer(oh, delimiter='\t')
    w.writerow(['ASV','domain','phylum','class','order','family','genus','species'])
    for asv in counts:
        r=pr2.get(asv,{})
        w.writerow([asv]+[clean_tax(r.get(k,'')) for k in ['domain','phylum','class','order','family','genus','species']])

with open(os.path.join(integrated_dir, 'taxonomy_JEDI_consensus.tsv'), 'w', newline='') as oh:
    w=csv.writer(oh, delimiter='\t')
    w.writerow(['ASV','Domain_JEDI','Phylum_JEDI','Class_JEDI','Order_JEDI','Family_JEDI','Genus_JEDI','Species_JEDI','Source'])
    for asv in counts:
        s=silva.get(asv,{k:'' for k in ['domain','phylum','class','order','family','genus','species']})
        p=pr2.get(asv,{k:'' for k in ['domain','phylum','class','order','family','genus','species']})
        domain=choose_domain(s,p)
        if domain in ('Bacteria','Archaea','Mitochondria','Plastid'):
            src='SILVA'
            base=s
        elif domain == 'Eukaryota':
            src='PR2'
            base=p
        else:
            src='UNRESOLVED'
            base=s if any(clean_tax(v) for v in s.values()) else p
        row=[asv, domain] + [clean_tax(base.get(k,'')) for k in ['phylum','class','order','family','genus','species']] + [src]
        consensus.append(row)
        w.writerow(row)
        total=sum(counts[asv].values())
        domain_summary[domain]+=total
        for sample,val in counts[asv].items():
            domain_per_sample[sample][domain]+=val

with open(os.path.join(integrated_dir, 'ASV_table_JEDI_counts.tsv'), 'w', newline='') as oh:
    w=csv.writer(oh, delimiter='\t')
    w.writerow(['ASV']+samples)
    for asv in counts:
        w.writerow([asv]+[counts[asv][s] for s in samples])

with open(os.path.join(integrated_dir, 'ASV_table_JEDI_with_taxonomy.tsv'), 'w', newline='') as oh:
    w=csv.writer(oh, delimiter='\t')
    w.writerow(['ASV','Domain_JEDI','Phylum_JEDI','Class_JEDI','Order_JEDI','Family_JEDI','Genus_JEDI','Species_JEDI','Source']+samples+['Sequence'])
    cdict={r[0]:r[1:] for r in consensus}
    for asv in counts:
        row=cdict[asv]
        w.writerow([asv]+row+[counts[asv][s] for s in samples]+[seqs.get(asv,'')])

with open(os.path.join(integrated_dir, 'domain_summary.tsv'), 'w', newline='') as oh:
    w=csv.writer(oh, delimiter='\t')
    w.writerow(['Domain','Total_reads'])
    for dom,total in sorted(domain_summary.items()):
        w.writerow([dom,total])

all_domains=sorted({d for smp in domain_per_sample.values() for d in smp.keys()})
with open(os.path.join(integrated_dir, 'domain_counts_per_sample.tsv'), 'w', newline='') as oh:
    w=csv.writer(oh, delimiter='\t')
    w.writerow(['Sample']+all_domains)
    for s in samples:
        w.writerow([s]+[domain_per_sample[s].get(d,0) for d in all_domains])

with open(os.path.join(integrated_dir, 'domain_relative_abundance_per_sample.tsv'), 'w', newline='') as oh:
    w=csv.writer(oh, delimiter='\t')
    w.writerow(['Sample']+all_domains)
    for s in samples:
        total=sum(domain_per_sample[s].values())
        vals=[(domain_per_sample[s].get(d,0)/total if total else 0) for d in all_domains]
        w.writerow([s]+vals)

for dom in all_domains:
    path=os.path.join(integrated_dir,'tables_by_domain',f'{dom}_ASV_table.tsv')
    dom_asvs=[r[0] for r in consensus if r[1]==dom]
    with open(path,'w',newline='') as oh:
        w=csv.writer(oh, delimiter='\t')
        w.writerow(['ASV']+samples)
        for asv in dom_asvs:
            w.writerow([asv]+[counts[asv][s] for s in samples])

with open(os.path.join(integrated_dir, 'alpha_diversity_unrarefied.tsv'), 'w', newline='') as oh:
    w=csv.writer(oh, delimiter='\t')
    w.writerow(['Sample','Observed_ASVs','Shannon','Simpson','Reads'])
    for s in samples:
        vec=[counts[a][s] for a in counts if counts[a][s] > 0]
        reads=sum(vec)
        observed=len(vec)
        if reads:
            ps=[v/reads for v in vec]
            sh=-sum(p*math.log(p) for p in ps if p>0)
            sim=1-sum(p*p for p in ps)
        else:
            sh=0.0
            sim=0.0
        w.writerow([s,observed,sh,sim,reads])

with open(os.path.join(integrated_dir, 'rarefaction_expected_ASVs.tsv'), 'w', newline='') as oh:
    w=csv.writer(oh, delimiter='\t')
    depths=[]
    totals={s:sum(counts[a][s] for a in counts) for s in samples}
    max_depth=max(totals.values()) if totals else 0
    if max_depth == 0:
        depths=[0]
    else:
        step=max(1000, max_depth//20)
        depths=list(range(step, max_depth+1, step))
    w.writerow(['Sample','Depth','Expected_ASVs'])
    for s in samples:
        total=totals[s]
        abund=[counts[a][s] for a in counts if counts[a][s] > 0]
        for m in depths:
            if total == 0 or m == 0:
                exp=0.0
            elif m >= total:
                exp=float(len(abund))
            else:
                exp=0.0
                for n in abund:
                    if total - n >= m:
                        num=math.lgamma(total-n+1)-math.lgamma(m+1)-math.lgamma(total-n-m+1)
                        den=math.lgamma(total+1)-math.lgamma(m+1)-math.lgamma(total-m+1)
                        p0=math.exp(num-den)
                    else:
                        p0=0.0
                    exp += 1-p0
            w.writerow([s,m,exp])
PY
}

write_manifest() {
  log "Écriture du manifest final"
  cat > "${RUN_MANIFEST}" <<EOF_MANIFEST
key	value
script_version	${SCRIPT_VERSION}
date	$(date '+%F %T')
project_root	${PROJECT_ROOT}
raw_dir	${RAW_DIR}
jedi_dir	${JEDI_DIR}
silva_dir	${SILVA_DIR}
pr2_dir	${PR2_DIR}
integrated_dir	${INTEGRATED_DIR}
nfcore_version	${NFCORE_VERSION}
profile	${PROFILE}
forward_primer	${FORWARD_PRIMER}
reverse_primer	${REVERSE_PRIMER}
trunc_len_f	${TRUNC_LEN_F}
trunc_len_r	${TRUNC_LEN_R}
mergepairs_strategy	consensus
silva_reference	${SILVA_REF}
pr2_reference	${PR2_REF}
preloaded_container	${PRELOAD_IMAGE_SIF}
nxf_singularity_cachedir	${NXF_SINGULARITY_CACHEDIR}
singularity_cachedir	${SINGULARITY_CACHEDIR}
singularity_tmpdir	${SINGULARITY_TMPDIR}
EOF_MANIFEST
}

validate_final_outputs() {
  local required=(
    "${SILVA_DIR}/dada2/ASV_seqs.fasta"
    "${SILVA_DIR}/dada2/table.tsv"
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
  local f
  for f in "${required[@]}"; do
    [[ -s "${f}" ]] || die "Sortie finale absente ou vide : ${f}"
  done
  printf 'SUCCESS\t%s\tJEDI pipeline completed\n' "$(date '+%F %T')" > "${SUCCESS_FLAG}"
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

  if [[ "${FORCE_SILVA}" == "1" || ! has_silva_success ]]; then
    log "ÉTAPE 1/3 : DADA2/JEDI et taxonomie SILVA"
    write_silva_params
    run_nfcore "${SILVA_PARAMS}" "${JEDI_DIR}/work/silva" "SILVA"
  else
    log "ÉTAPE 1/3 : SILVA déjà disponible, réutilisation"
  fi

  has_silva_success || die "Les sorties SILVA minimales ne sont pas présentes après exécution"

  if [[ "${FORCE_PR2}" == "1" || ! has_pr2_success ]]; then
    log "ÉTAPE 2/3 : taxonomie PR2 sur ASV existants"
    rm -rf "${PR2_DIR}" && mkdir -p "${PR2_DIR}"
    write_pr2_params
    run_nfcore "${PR2_PARAMS}" "${JEDI_DIR}/work/pr2" "PR2"
  else
    log "ÉTAPE 2/3 : PR2 déjà disponible, réutilisation"
  fi

  has_pr2_success || die "Les sorties PR2 minimales ne sont pas présentes après exécution"

  log "ÉTAPE 3/3 : intégration cross-domain"
  integrate_outputs
  write_manifest
  validate_final_outputs
  log "Pipeline JEDI terminé avec succès"
}

main "$@"
