#!/usr/bin/env bash
# ==============================================================================
# Pipeline JEDI complet et independant pour le projet valormicro_V2
#
# Principe :
#   FASTQ bruts -> Cutadapt -> DADA2 -> consensus merge pairs -> ASV
#   -> taxonomie SILVA (Bacteria/Archaea/plastides)
#   -> taxonomie PR2 (Eukaryota)
#   -> integration cross-domain Bacteria/Archaea/Eukaryota/organites
#
# IMPORTANT : ce pipeline est adapte aux amplicons JEDI 515F-Y/926R.
# Il ne peut pas creer de signal eucaryote/archeen absent de la PCR initiale.
#
# Lancement conseille :
#   nohup bash jedi_valormicro.sh > jedi_launcher.log 2>&1 &
#
# Reprise apres interruption : relancer exactement la meme commande ; Nextflow
# reutilisera son cache avec -resume.
# ==============================================================================

set -Eeuo pipefail
shopt -s nullglob
IFS=$'\n\t'

# ------------------------------- CONFIGURATION -------------------------------
PROJECT_DIR="${PROJECT_DIR:-/nvme/bio/data_fungi/valormicro_V2}"
RAW_DIR="${RAW_DIR:-${PROJECT_DIR}/01_raw_data}"
INFO_XLSX="${INFO_XLSX:-${RAW_DIR}/00_infos_data.xlsx}"

# Dossier volontairement distinct de 02_amplicon_pipeline.
JEDI_ROOT="${JEDI_ROOT:-${PROJECT_DIR}/03_JEDI_pipeline}"
INPUT_DIR="${JEDI_ROOT}/00_inputs"
SILVA_OUT="${JEDI_ROOT}/01_nfcore_silva"
PR2_OUT="${JEDI_ROOT}/02_nfcore_pr2"
INTEGRATED_OUT="${JEDI_ROOT}/03_integrated"
WORK_DIR="${JEDI_ROOT}/work"
LOG_DIR="${JEDI_ROOT}/logs"
DB_CACHE="${JEDI_ROOT}/reference_databases"
TMP_DIR="${JEDI_ROOT}/tmp"

# Version reproductible de nf-core/ampliseq. La strategie consensus est
# disponible depuis la version 2.13.
NFCORE_VERSION="${NFCORE_VERSION:-2.18.0}"

# Bases du papier JEDI, avec versions explicites et modifiables.
SILVA_DB="${SILVA_DB:-silva=138.2}"
PR2_DB="${PR2_DB:-pr2=5.0.0}"

# Amorces JEDI 515F-Y / 926R, ecrites 5' -> 3'.
FW_PRIMER="${FW_PRIMER:-GTGYCAGCMGCCGCGGTAA}"
RV_PRIMER="${RV_PRIMER:-CCGYCAATTYMTTTRAGTTT}"

# nf-core determine automatiquement les longueurs de troncature juste avant
# la chute de qualite mediane sous Q25. Pour imposer les longueurs du papier
# (par exemple 231/230 pour des reads 2x250), definir :
#   export JEDI_TRUNCLEN_F=231 JEDI_TRUNCLEN_R=230
#JEDI_TRUNCLEN_F="${JEDI_TRUNCLEN_F:-auto}"
#JEDI_TRUNCLEN_R="${JEDI_TRUNCLEN_R:-auto}"
export JEDI_TRUNCLEN_F=231
export JEDI_TRUNCLEN_R=230
TRUNC_QMIN="${TRUNC_QMIN:-25}"
TRUNC_RMIN="${TRUNC_RMIN:-0.75}"
MAX_EE="${MAX_EE:-2}"

# Ce seuil ne filtre que les analyses de diversite QIIME2, pas les tables ASV
# brutes. A ajuster apres inspection des courbes de rarefaction.
DIVERSITY_RAREFACTION_DEPTH="${DIVERSITY_RAREFACTION_DEPTH:-500}"

# Environnement utilise uniquement pour lire le tableur Excel.
EXCEL_ENV="${EXCEL_ENV:-excel_tools}"

# Profil d'execution : si NXF_PROFILE est vide, detection automatique.
NXF_PROFILE="${NXF_PROFILE:-}"

JEDI=/nvme/bio/data_fungi/valormicro_V2/03_JEDI_pipeline

mkdir -p "$JEDI/tmp/qiime2"
mkdir -p "$JEDI/tmp/singularity"
mkdir -p "$JEDI/tmp/xdg_config"
mkdir -p "$JEDI/tmp/mpl"
mkdir -p "$JEDI/tmp/numba"

chmod -R u+rwx "$JEDI/tmp"

# ------------------------------- FONCTIONS -----------------------------------
log() {
    printf '[%(%F %T)T] %s\n' -1 "$*"
}

die() {
    log "ERREUR : $*" >&2
    exit 1
}

on_error() {
    local rc=$?
    log "ERREUR : code ${rc}, ligne ${BASH_LINENO[0]} : ${BASH_COMMAND}" >&2
    exit "${rc}"
}
trap on_error ERR

choose_profile() {
    if [[ -n "${NXF_PROFILE}" ]]; then
        printf '%s' "${NXF_PROFILE}"
    elif command -v apptainer >/dev/null 2>&1; then
        printf '%s' "apptainer"
    elif command -v singularity >/dev/null 2>&1; then
        printf '%s' "singularity"
    elif command -v docker >/dev/null 2>&1; then
        printf '%s' "docker"
    elif command -v conda >/dev/null 2>&1; then
        printf '%s' "conda"
    else
        return 1
    fi
}

run_nextflow() {
    local params_file="$1"
    local work_subdir="$2"
    log "Commande : nextflow run nf-core/ampliseq -r ${NFCORE_VERSION} -profile ${PROFILE} -params-file ${params_file} -work-dir ${work_subdir} -resume"
   nextflow run nf-core/ampliseq \
    -r "${NFCORE_VERSION}" \
    -profile "${PROFILE}" \
    -params-file "${params_file}" \
    -work-dir "${work_subdir}" \
    -resume \
    -c "${JEDI_ROOT}/singularity_timeout.config"
}

cat > /nvme/bio/data_fungi/valormicro_V2/03_JEDI_pipeline/singularity_timeout.config <<'EOF'
singularity {
    pullTimeout = '2h'
}
EOF

# ----------------------------- PREPARATION -----------------------------------
mkdir -p "${INPUT_DIR}" "${SILVA_OUT}" "${PR2_OUT}" "${INTEGRATED_OUT}" \
         "${WORK_DIR}" "${LOG_DIR}" "${DB_CACHE}" "${TMP_DIR}"

exec > >(tee -a "${LOG_DIR}/jedi_pipeline.log") 2>&1

export TMPDIR="${TMP_DIR}"
export NXF_HOME="${JEDI_ROOT}/.nextflow"
export NXF_OPTS="${NXF_OPTS:--Xms1g -Xmx4g}"
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-${JEDI_ROOT}/container_cache}"
export NXF_SINGULARITY_PULL_TIMEOUT="2h"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${JEDI_ROOT}/container_cache}"
mkdir -p "${NXF_HOME}" "${NXF_SINGULARITY_CACHEDIR}"

command -v nextflow >/dev/null 2>&1 || die "Nextflow est introuvable dans le PATH."
[[ -f "${INFO_XLSX}" ]] || die "Tableur absent : ${INFO_XLSX}"
[[ -d "${RAW_DIR}" ]] || die "Dossier FASTQ absent : ${RAW_DIR}"

# Espaces temporaires persistants et inscriptibles pour QIIME2/Singularity.
export TMPDIR="/nvme/bio/data_fungi/valormicro_V2/03_JEDI_pipeline/tmp"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

export XDG_CONFIG_HOME="$TMPDIR/xdg_config"
export MPLCONFIGDIR="$TMPDIR/mpl"
export NUMBA_CACHE_DIR="$TMPDIR/numba"

mkdir -p "$TMPDIR/qiime2" "$XDG_CONFIG_HOME" "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"

cat > /nvme/bio/data_fungi/valormicro_V2/03_JEDI_pipeline/singularity_timeout.config <<'EOF'
singularity {
    enabled = true
    autoMounts = true
    pullTimeout = '2h'
    cacheDir = '/nvme/bio/data_fungi/valormicro_V2/03_JEDI_pipeline/container_cache'
    runOptions = '--bind /nvme/bio/data_fungi/valormicro_V2/03_JEDI_pipeline/tmp:/nvme/bio/data_fungi/valormicro_V2/03_JEDI_pipeline/tmp'
}

process {
    withName: 'NFCORE_AMPLISEQ:AMPLISEQ:FORMATTAXONOMY' {
        container = '/nvme/bio/data_fungi/valormicro_V2/03_JEDI_pipeline/container_cache/biocontainers-v1.2.0_cv1.sif'
    }
}
EOF

PROFILE="$(choose_profile)" || die "Aucun moteur disponible (Apptainer, Singularity, Docker ou Conda)."
log "Dossier JEDI independant : ${JEDI_ROOT}"
log "nf-core/ampliseq : ${NFCORE_VERSION}"
log "Profil : ${PROFILE}"
log "Bases : ${SILVA_DB} et ${PR2_DB}"
log "Amorces : ${FW_PRIMER} / ${RV_PRIMER}"
nextflow -version

# Python avec pandas/openpyxl : environnement historique si disponible, sinon
# Chargement explicite de Conda, necessaire pour l'environnement excel_tools.
command -v conda >/dev/null 2>&1 \
    || die "Conda est introuvable dans le PATH."

CONDA_BASE="$(conda info --base 2>/dev/null || true)"
[[ -n "${CONDA_BASE}" && -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]] \
    || die "Impossible de localiser conda.sh."

source "${CONDA_BASE}/etc/profile.d/conda.sh"

conda activate "${EXCEL_ENV}" \
    || die "Impossible d'activer l'environnement ${EXCEL_ENV}."

command -v python >/dev/null 2>&1 \
    || die "Python absent de l'environnement ${EXCEL_ENV}."

python -c 'import pandas, openpyxl' \
    || die "pandas et/ou openpyxl absents de l'environnement ${EXCEL_ENV}."

PYTHON_CMD=(python)

SAMPLESHEET="${INPUT_DIR}/samplesheet_jedi.tsv"
METADATA="${INPUT_DIR}/metadata_jedi.tsv"
ID_MAP="${INPUT_DIR}/sample_id_mapping.tsv"
METADATA_FLAG="${INPUT_DIR}/metadata_is_informative.flag"
CONTROL_FLAG="${INPUT_DIR}/controls_detected.flag"

log "Construction des entrees nf-core depuis ${INFO_XLSX}"
"${PYTHON_CMD[@]}" - "${INFO_XLSX}" "${RAW_DIR}" "${SAMPLESHEET}" \
"${METADATA}" "${ID_MAP}" "${METADATA_FLAG}" "${CONTROL_FLAG}" <<'PY'
import re
import sys
import unicodedata
from pathlib import Path

import pandas as pd

xlsx, raw_dir, samplesheet, metadata, id_map, metadata_flag, control_flag = sys.argv[1:]
raw_dir = Path(raw_dir).resolve()
df = pd.read_excel(xlsx, dtype=str)
df.columns = [str(c).strip() for c in df.columns]
required = ["R1", "R2", "New_label"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit(f"Colonnes obligatoires absentes : {', '.join(missing)}")

for c in df.columns:
    df[c] = df[c].fillna("").astype(str).str.strip()

any_required = (df["R1"] != "") | (df["R2"] != "") | (df["New_label"] != "")
partial = df.loc[any_required & ((df["R1"] == "") | (df["R2"] == "") | (df["New_label"] == "")), required]
if not partial.empty:
    raise SystemExit("Lignes incompletes R1/R2/New_label :\n" + partial.to_string(index=False))
df = df.loc[(df["R1"] != "") & (df["R2"] != "") & (df["New_label"] != "")].copy()
if df.empty:
    raise SystemExit("Aucun echantillon complet dans le tableur.")

def ascii_text(value):
    return unicodedata.normalize("NFKD", str(value)).encode("ascii", "ignore").decode()

def sample_id(value):
    value = re.sub(r"[^A-Za-z0-9_]+", "_", ascii_text(value).strip())
    value = re.sub(r"_+", "_", value).strip("_")
    if not value:
        raise SystemExit("Un New_label ne donne aucun identifiant valide.")
    if not value[0].isalpha():
        value = "S_" + value
    return value[:36]

df["sample"] = [sample_id(x) for x in df["New_label"]]
if df["sample"].duplicated().any():
    dup = df.loc[df["sample"].duplicated(keep=False), ["New_label", "sample"]]
    raise SystemExit("IDs dupliques apres normalisation :\n" + dup.to_string(index=False))

def fq_path(value):
    p = Path(value)
    if not p.is_absolute():
        p = raw_dir / p
    p = p.resolve()
    if not p.is_file() or p.stat().st_size == 0:
        raise SystemExit(f"FASTQ absent ou vide : {p}")
    if not (str(p).endswith(".fastq.gz") or str(p).endswith(".fq.gz")):
        raise SystemExit(f"FASTQ non compresse ou extension invalide : {p}")
    return str(p)

out = pd.DataFrame({
    "sample": df["sample"],
    "fastq_1": [fq_path(x) for x in df["R1"]],
    "fastq_2": [fq_path(x) for x in df["R2"]],
})

# Colonnes optionnelles reconnues uniquement par nom explicite, sans inventer
# le statut des controles.
norm_to_original = {re.sub(r"[^a-z0-9]+", "", ascii_text(c).lower()): c for c in df.columns}
run_col = next((norm_to_original[k] for k in ("run", "sequencingrun", "runid") if k in norm_to_original), None)
control_col = next((norm_to_original[k] for k in ("control", "controle", "negativecontrol") if k in norm_to_original), None)
quant_col = next((norm_to_original[k] for k in ("quantreading", "dnaquantity", "dnaconcentration") if k in norm_to_original), None)

if run_col and (df[run_col] != "").any():
    out["run"] = df[run_col].replace("", "run1")

if control_col:
    def control_value(x):
        y = ascii_text(x).strip().lower()
        controls = {"control", "controle", "negative", "negatif", "blank", "blanc", "ntc"}
        samples = {"sample", "echantillon", "positive", "positif"}
        if y in controls or any(t in y for t in ("negative control", "controle negatif", "blank")):
            return "control"
        if y in samples or y == "":
            return "sample"
        raise SystemExit(f"Valeur de controle non interpretable : {x!r}")
    out["control"] = [control_value(x) for x in df[control_col]]
    if (out["control"] == "control").any():
        Path(control_flag).write_text("controls present\n")

if quant_col and (df[quant_col] != "").all():
    quant = pd.to_numeric(df[quant_col].str.replace(",", ".", regex=False), errors="coerce")
    if quant.notna().all():
        out["quant_reading"] = quant

out.to_csv(samplesheet, sep="\t", index=False)
pd.DataFrame({"New_label_original": df["New_label"], "sample": df["sample"]}).to_csv(id_map, sep="\t", index=False)

# Metadonnees : retrait des chemins techniques et normalisation des entetes.
meta = df.drop(columns=["R1", "R2", "New_label", "sample"], errors="ignore").copy()

def column_name(value):
    value = re.sub(r"[^A-Za-z0-9_]+", "_", ascii_text(value).strip())
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "metadata"

meta.columns = [column_name(c) for c in meta.columns]
meta.insert(0, "ID", df["sample"].values)
meta.to_csv(metadata, sep="\t", index=False)

# nf-core utilise les metadonnees pour les analyses aval seulement si au moins
# une variable categorielle contient plusieurs groupes et n'est pas unique par
# echantillon.
informative = False
for c in meta.columns[1:]:
    vals = meta[c].astype(str)
    n = vals[vals != ""].nunique()
    if 2 <= n < len(meta):
        informative = True
        break
if informative:
    Path(metadata_flag).write_text("informative metadata present\n")

print(f"{len(out)} echantillons ecrits dans {samplesheet}")
print(f"Metadonnees informatives : {informative}")
print(f"Controles detectes : {Path(control_flag).exists()}")
PY

[[ -s "${SAMPLESHEET}" ]] || die "Samplesheet non produit."
[[ -s "${METADATA}" ]] || die "Metadata non produit."

log "Verification gzip de tous les FASTQ declares"
while IFS=$'\t' read -r sample fastq_1 fastq_2 rest; do
    [[ "${sample}" == "sample" ]] && continue
    gzip -t "${fastq_1}"
    gzip -t "${fastq_2}"
done < "${SAMPLESHEET}"

# Copie du script lance pour provenance.
cp -f "$0" "${INPUT_DIR}/pipeline_jedi_executed.sh" 2>/dev/null || true

# ------------------------ PARAMETRES NF-CORE PRINCIPAUX ----------------------
SILVA_PARAMS="${INPUT_DIR}/params_silva.yaml"
cat > "${SILVA_PARAMS}" <<YAML
input: "${SAMPLESHEET}"
outdir: "${SILVA_OUT}"
FW_primer: "${FW_PRIMER}"
RV_primer: "${RV_PRIMER}"
ref_taxonomy_storage: "${DB_CACHE}"
save_intermediates: true

# JEDI : fusion des paires chevauchantes et concatenation des non-chevauchantes.
mergepairs_strategy: "consensus"
mergepairs_consensus_match: 1
mergepairs_consensus_mismatch: -2
mergepairs_consensus_gap: -4
mergepairs_consensus_minoverlap: 12
mergepairs_consensus_maxmismatch: 0
mergepairs_consensus_percentile_cutoff: 0.001

# Qualite DADA2.
trunc_qmin: ${TRUNC_QMIN}
trunc_rmin: ${TRUNC_RMIN}
max_ee: ${MAX_EE}
min_len: 50
sample_inference: "independent"

# Taxonomie prokaryote et organites.
dada_ref_taxonomy: "${SILVA_DB}"
cut_dada_ref_taxonomy: true
skip_dada_addspecies: true

# Ne perdre aucun domaine ni ASV rare dans les tables principales.
exclude_taxa: "none"
min_frequency: 1
min_samples: 1
diversity_rarefaction_depth: ${DIVERSITY_RAREFACTION_DEPTH}
report_title: "JEDI valormicro - SILVA, PR2 et consensus merge pairs"
YAML

if [[ -f "${METADATA_FLAG}" ]]; then
    printf 'metadata: "%s"\n' "${METADATA}" >> "${SILVA_PARAMS}"
else
    log "Aucune variable categorielle informative : analyses dependantes des metadonnees ignorees."
fi

if [[ "${JEDI_TRUNCLEN_F}" != "auto" || "${JEDI_TRUNCLEN_R}" != "auto" ]]; then
    [[ "${JEDI_TRUNCLEN_F}" =~ ^[0-9]+$ && "${JEDI_TRUNCLEN_R}" =~ ^[0-9]+$ ]] \
        || die "JEDI_TRUNCLEN_F et JEDI_TRUNCLEN_R doivent etre deux entiers ou tous deux 'auto'."
    printf 'trunclenf: %s\ntrunclenr: %s\n' "${JEDI_TRUNCLEN_F}" "${JEDI_TRUNCLEN_R}" >> "${SILVA_PARAMS}"
fi

if [[ -f "${CONTROL_FLAG}" ]]; then
    log "Controles negatifs detectes : decontam sera calcule mais pas applique aux tables finales."
fi

# ------------------------- EXECUTION JEDI + SILVA ----------------------------
log "ETAPE 1/3 : inference des ASV JEDI et classification SILVA"
run_nextflow "${SILVA_PARAMS}" "${WORK_DIR}/silva"

ASV_FASTA="${SILVA_OUT}/dada2/ASV_seqs.fasta"
ASV_TABLE="${SILVA_OUT}/dada2/ASV_table.tsv"
[[ -s "${ASV_FASTA}" ]] || die "ASV FASTA absent apres nf-core : ${ASV_FASTA}"
[[ -s "${ASV_TABLE}" ]] || die "Table ASV absente apres nf-core : ${ASV_TABLE}"

# ---------------------------- TAXONOMIE PR2 ----------------------------------
PR2_PARAMS="${INPUT_DIR}/params_pr2.yaml"
cat > "${PR2_PARAMS}" <<YAML
input_fasta: "${ASV_FASTA}"
outdir: "${PR2_OUT}"
FW_primer: "${FW_PRIMER}"
RV_primer: "${RV_PRIMER}"
ref_taxonomy_storage: "${DB_CACHE}"
save_intermediates: true

dada_ref_taxonomy: "${PR2_DB}"
cut_dada_ref_taxonomy: true
skip_dada_addspecies: true
skip_qiime: true
report_title: "JEDI valormicro - classification eucaryote PR2"
YAML

log "ETAPE 2/3 : classification des memes ASV avec PR2"
run_nextflow "${PR2_PARAMS}" "${WORK_DIR}/pr2"

# ---------------------- INTEGRATION CROSS-DOMAIN -----------------------------
log "ETAPE 3/3 : integration SILVA + PR2 et tables par domaine"
"${PYTHON_CMD[@]}" - "${SILVA_OUT}" "${PR2_OUT}" "${ASV_TABLE}" \
    "${ASV_FASTA}" "${INTEGRATED_OUT}" <<'PY'
import csv
import glob
import shutil
import sys
from collections import defaultdict
from pathlib import Path

silva_out, pr2_out, table_path, fasta_path, outdir = map(Path, sys.argv[1:])
outdir.mkdir(parents=True, exist_ok=True)

def taxonomy_file(root, label):
    candidates = []
    for pattern in ("dada2/ASV_tax.*.tsv", "dada2/ASV_tax.tsv"):
        for p in root.glob(pattern):
            name = p.name.lower()
            if "species" not in name and p.is_file() and p.stat().st_size > 0:
                candidates.append(p)
    if not candidates:
        raise SystemExit(f"Taxonomie {label} introuvable dans {root}/dada2")
    candidates = sorted(set(candidates), key=lambda p: (label.lower() not in p.name.lower(), len(p.name)))
    return candidates[0]

def read_taxonomy(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"Taxonomie sans entete : {path}")
        id_col = next((c for c in reader.fieldnames if c.lower() in {"asv_id", "feature id", "featureid", "id"}), reader.fieldnames[0])
        rank_cols = [c for c in reader.fieldnames if c != id_col and c.lower() not in {"sequence", "confidence"}]
        data = {}
        for row in reader:
            asv = row[id_col]
            ranks = [(row.get(c) or "").strip() for c in rank_cols]
            data[asv] = (rank_cols, ranks)
    return data

def read_fasta(path):
    seqs, current = {}, None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current:
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}

def clean(x):
    return x.strip().strip(";")

def assigned(ranks):
    bad = {"", "na", "nan", "unassigned", "unknown", "root"}
    return [clean(x) for x in ranks if clean(x).lower() not in bad]

def joined(ranks):
    return ";".join(assigned(ranks))

def has(text, words):
    low = text.lower()
    return any(w in low for w in words)

silva_path = taxonomy_file(silva_out, "silva")
pr2_path = taxonomy_file(pr2_out, "pr2")
silva = read_taxonomy(silva_path)
pr2 = read_taxonomy(pr2_path)
seqs = read_fasta(fasta_path)

with table_path.open(encoding="utf-8-sig", newline="") as handle:
    rows = list(csv.reader(handle, delimiter="\t"))
if not rows:
    raise SystemExit("Table ASV vide")
header = rows[0]
id_index = 0
sample_names = header[1:]

consensus = {}
for row in rows[1:]:
    if not row:
        continue
    asv = row[id_index]
    s_ranks = silva.get(asv, ([], []))[1]
    p_ranks = pr2.get(asv, ([], []))[1]
    s = joined(s_ranks)
    p = joined(p_ranks)
    sl, pl = s.lower(), p.lower()

    # Les organites sont testes avant Bacteria car les plastides sont places
    # dans les Cyanobacteria par les references 16S.
    if has(sl, ["chloroplast", "plastid"]):
        domain, source, tax = "Eukaryota_plastid", "SILVA", s
    elif "mitochond" in sl:
        domain, source, tax = "Eukaryota_mitochondria", "SILVA", s
    elif has(sl, ["archaea", "d_0__archaea", "k__archaea"]):
        domain, source, tax = "Archaea", "SILVA", s
    elif has(sl, ["bacteria", "d_0__bacteria", "k__bacteria"]):
        domain, source, tax = "Bacteria", "SILVA", s
    elif has(pl, ["eukaryota", "eukaryote"]):
        domain, source, tax = "Eukaryota", "PR2", p
    elif p:
        # PR2 est une base eucaryote ; une annotation informative PR2 est
        # conservee mais explicitement marquee pour audit.
        domain, source, tax = "Eukaryota_candidate", "PR2", p
    elif s:
        domain, source, tax = "Unresolved_SILVA", "SILVA", s
    else:
        domain, source, tax = "Unassigned", "none", "Unassigned"

    consensus[asv] = {
        "domain": domain,
        "source": source,
        "taxonomy": tax or "Unassigned",
        "silva": s or "Unassigned",
        "pr2": p or "Unassigned",
        "length": len(seqs.get(asv, "")),
    }

# Taxonomie auditable, les deux annotations restant toujours visibles.
tax_out = outdir / "taxonomy_JEDI_consensus.tsv"
with tax_out.open("w", encoding="utf-8", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["ASV_ID", "JEDI_domain", "selected_source", "selected_taxonomy", "SILVA_taxonomy", "PR2_taxonomy", "ASV_length"])
    for asv in sorted(consensus):
        x = consensus[asv]
        w.writerow([asv, x["domain"], x["source"], x["taxonomy"], x["silva"], x["pr2"], x["length"]])

# Table brute conservee, enrichie sans rarefaction ni filtrage de prevalence.
joined_out = outdir / "ASV_table_JEDI_with_taxonomy.tsv"
counts_out = outdir / "ASV_table_JEDI_counts.tsv"
with joined_out.open("w", encoding="utf-8", newline="") as full, counts_out.open("w", encoding="utf-8", newline="") as counts:
    wf = csv.writer(full, delimiter="\t", lineterminator="\n")
    wc = csv.writer(counts, delimiter="\t", lineterminator="\n")
    wc.writerow(header)
    wf.writerow(["ASV_ID", "JEDI_domain", "selected_source", "selected_taxonomy", "SILVA_taxonomy", "PR2_taxonomy", "ASV_length"] + sample_names)
    for row in rows[1:]:
        if not row:
            continue
        asv = row[0]
        x = consensus[asv]
        wc.writerow(row)
        wf.writerow([asv, x["domain"], x["source"], x["taxonomy"], x["silva"], x["pr2"], x["length"]] + row[1:])

# Comptages et abondances relatives par domaine et par echantillon.
domain_counts = defaultdict(lambda: [0] * len(sample_names))
domain_asvs = defaultdict(int)
for row in rows[1:]:
    if not row:
        continue
    asv = row[0]
    domain = consensus[asv]["domain"]
    values = [int(float(x or 0)) for x in row[1:]]
    domain_asvs[domain] += 1
    domain_counts[domain] = [a + b for a, b in zip(domain_counts[domain], values)]

domains = sorted(domain_counts)
with (outdir / "domain_counts_per_sample.tsv").open("w", encoding="utf-8", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["JEDI_domain"] + sample_names)
    for d in domains:
        w.writerow([d] + domain_counts[d])

sample_totals = [sum(domain_counts[d][i] for d in domains) for i in range(len(sample_names))]
with (outdir / "domain_relative_abundance_per_sample.tsv").open("w", encoding="utf-8", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["JEDI_domain"] + sample_names)
    for d in domains:
        vals = [(domain_counts[d][i] / sample_totals[i]) if sample_totals[i] else 0 for i in range(len(sample_names))]
        w.writerow([d] + [f"{x:.10f}" for x in vals])

with (outdir / "domain_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["JEDI_domain", "ASV_richness", "total_reads"])
    for d in domains:
        w.writerow([d, domain_asvs[d], sum(domain_counts[d])])

shutil.copy2(fasta_path, outdir / "ASV_sequences_JEDI.fasta")
shutil.copy2(silva_path, outdir / "taxonomy_SILVA_original.tsv")
shutil.copy2(pr2_path, outdir / "taxonomy_PR2_original.tsv")

print(f"SILVA : {silva_path}")
print(f"PR2 : {pr2_path}")
print(f"ASV integres : {len(consensus)}")
print(f"Domaines : {', '.join(domains)}")
PY

# Provenance finale.
{
    printf 'date\t%s\n' "$(date --iso-8601=seconds)"
    printf 'project_dir\t%s\n' "${PROJECT_DIR}"
    printf 'jedi_root\t%s\n' "${JEDI_ROOT}"
    printf 'nfcore_ampliseq\t%s\n' "${NFCORE_VERSION}"
    printf 'profile\t%s\n' "${PROFILE}"
    printf 'forward_primer\t%s\n' "${FW_PRIMER}"
    printf 'reverse_primer\t%s\n' "${RV_PRIMER}"
    printf 'silva_database\t%s\n' "${SILVA_DB}"
    printf 'pr2_database\t%s\n' "${PR2_DB}"
    printf 'mergepairs_strategy\tconsensus\n'
    printf 'trunclen_f\t%s\n' "${JEDI_TRUNCLEN_F}"
    printf 'trunclen_r\t%s\n' "${JEDI_TRUNCLEN_R}"
    printf 'trunc_qmin\t%s\n' "${TRUNC_QMIN}"
    printf 'max_ee\t%s\n' "${MAX_EE}"
} > "${INTEGRATED_OUT}/run_manifest.tsv"

log "Pipeline JEDI termine avec succes."
log "Rapport principal : ${SILVA_OUT}/summary_report/summary_report.html"
log "MultiQC : ${SILVA_OUT}/multiqc/multiqc_report.html"
log "ASV + taxonomie integree : ${INTEGRATED_OUT}/ASV_table_JEDI_with_taxonomy.tsv"
log "Resume par domaine : ${INTEGRATED_OUT}/domain_summary.tsv"
log "Courbes de rarefaction : ${SILVA_OUT}/qiime2/alpha-rarefaction/index.html"
