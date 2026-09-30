#!/usr/bin/env bash
# ==============================================================================
# Pipeline JEDI valormicro_V2 — version corrigée et relançable
#
# Objectif :
#   FASTQ paired-end -> Cutadapt -> DADA2/JEDI consensus -> ASV
#   -> SILVA 138.2 (Bacteria, Archaea, organites)
#   -> PR2 5.0.0 (Eukaryota)
#   -> intégration cross-domain, tables par domaine, alpha-diversité et
#      raréfaction analytique.
#
# Cette version désactive les analyses aval QIIME2/Emperor dans nf-core : elles
# ne sont pas nécessaires à l'inférence/classification et sont la cause exacte
# de l'arrêt observé dans 50_jedi_valormicro-2.out. Les tables aval sont créées
# ici de manière déterministe après les deux classifications.
#
# Relance : la commande est identique ; les résultats valides et le cache
# Nextflow sont réutilisés. Pour forcer une étape :
#   FORCE_SILVA=1 bash 50_jedi_valormicro_fixed.sh
#   FORCE_PR2=1   bash 50_jedi_valormicro_fixed.sh
#   FORCE_ALL=1   bash 50_jedi_valormicro_fixed.sh
# ==============================================================================

set -Eeuo pipefail
shopt -s nullglob
IFS=$'\n\t'
umask 002

# -------------------------------- CONFIGURATION ------------------------------
PROJECT_DIR="${PROJECT_DIR:-/nvme/bio/data_fungi/valormicro_V2}"
RAW_DIR="${RAW_DIR:-${PROJECT_DIR}/01_raw_data}"
INFO_XLSX="${INFO_XLSX:-${RAW_DIR}/00_infos_data.xlsx}"
JEDI_ROOT="${JEDI_ROOT:-${PROJECT_DIR}/03_JEDI_pipeline}"

INPUT_DIR="${JEDI_ROOT}/00_inputs"
SILVA_OUT="${JEDI_ROOT}/01_nfcore_silva"
PR2_OUT="${JEDI_ROOT}/02_nfcore_pr2"
INTEGRATED_OUT="${JEDI_ROOT}/03_integrated"
WORK_DIR="${JEDI_ROOT}/work"
LOG_DIR="${JEDI_ROOT}/logs"
DB_CACHE="${JEDI_ROOT}/reference_databases"
CONTAINER_CACHE="${JEDI_ROOT}/container_cache"
TMP_ROOT="${JEDI_ROOT}/tmp"
LAUNCH_DIR="${CONTAINER_CACHE}"

NFCORE_VERSION="${NFCORE_VERSION:-2.18.0}"
SILVA_DB="${SILVA_DB:-silva=138.2}"
PR2_DB="${PR2_DB:-pr2=5.0.0}"

# Amorces JEDI 515F-Y / 926R (5' -> 3').
FW_PRIMER="${FW_PRIMER:-GTGYCAGCMGCCGCGGTAA}"
RV_PRIMER="${RV_PRIMER:-CCGYCAATTYMTTTRAGTTT}"

# Valeurs utilisées dans le protocole JEDI sur des reads 2 x 250.
# Mettre 0/0 pour ne pas tronquer, ou modifier via variables d'environnement.
TRUNCLEN_F="${TRUNCLEN_F:-231}"
TRUNCLEN_R="${TRUNCLEN_R:-230}"
TRUNC_QMIN="${TRUNC_QMIN:-25}"
TRUNC_RMIN="${TRUNC_RMIN:-0.75}"
MAX_EE="${MAX_EE:-2}"
MIN_LEN="${MIN_LEN:-50}"

# Nombre de profondeurs de raréfaction analytique par échantillon.
RAREFACTION_POINTS="${RAREFACTION_POINTS:-20}"

EXCEL_ENV="${EXCEL_ENV:-excel_tools}"
NXF_PROFILE="${NXF_PROFILE:-}"
FORCE_SILVA="${FORCE_SILVA:-0}"
FORCE_PR2="${FORCE_PR2:-0}"
FORCE_ALL="${FORCE_ALL:-0}"

# --------------------------------- FONCTIONS ---------------------------------
log() {
    printf '[%(%F %T)T] %s\n' -1 "$*"
}

die() {
    log "ERREUR : $*" >&2
    exit 1
}

on_error() {
    local rc=$?
    local line="${BASH_LINENO[0]:-inconnue}"
    log "ERREUR : code ${rc}, ligne ${line}, commande : ${BASH_COMMAND}" >&2
    exit "$rc"
}
trap on_error ERR

choose_profile() {
    if [[ -n "$NXF_PROFILE" ]]; then
        printf '%s' "$NXF_PROFILE"
    elif command -v apptainer >/dev/null 2>&1; then
        printf '%s' apptainer
    elif command -v singularity >/dev/null 2>&1; then
        printf '%s' singularity
    elif command -v docker >/dev/null 2>&1; then
        printf '%s' docker
    else
        return 1
    fi
}

is_uint() {
    [[ "$1" =~ ^[0-9]+$ ]]
}

has_taxonomy() {
    local root="$1"
    local found
    found="$(find "${root}/dada2" -maxdepth 1 -type f \
        -name 'ASV_tax.*.tsv' ! -name '*species*' -size +0c \
        -print -quit 2>/dev/null || true)"
    [[ -n "$found" ]]
}

silva_core_complete() {
    [[ -s "${SILVA_OUT}/dada2/ASV_seqs.fasta" ]] &&
    [[ -s "${SILVA_OUT}/dada2/ASV_table.tsv" ]] &&
    has_taxonomy "$SILVA_OUT"
}

pr2_core_complete() {
    has_taxonomy "$PR2_OUT"
}

run_nextflow() {
    local params_file="$1"
    local work_subdir="$2"
    local run_label="$3"

    log "Lancement nf-core/ampliseq (${run_label})"
    log "Paramètres : ${params_file}"
    log "Work      : ${work_subdir}"

    (
        cd "$LAUNCH_DIR"
        nextflow run nf-core/ampliseq \
            -r "$NFCORE_VERSION" \
            -profile "$PROFILE" \
            -params-file "$params_file" \
            -work-dir "$work_subdir" \
            -c "${INPUT_DIR}/nextflow_infrastructure.config" \
            -resume \
            -ansi-log false
    )
}

# ------------------------------- PRÉPARATION ---------------------------------
mkdir -p "$INPUT_DIR" "$SILVA_OUT" "$PR2_OUT" "$INTEGRATED_OUT" \
    "$WORK_DIR" "$LOG_DIR" "$DB_CACHE" "$CONTAINER_CACHE" \
    "$TMP_ROOT/general" "$TMP_ROOT/nextflow" "$TMP_ROOT/xdg" \
    "$TMP_ROOT/mpl" "$TMP_ROOT/numba"

exec > >(tee -a "${LOG_DIR}/jedi_pipeline_fixed.log") 2>&1

# Verrou simple contre deux lancements simultanés.
LOCK_DIR="${JEDI_ROOT}/.jedi_fixed.lock"
if ! mkdir "$LOCK_DIR" 2>/dev/null; then
    die "Un autre lancement semble actif (${LOCK_DIR}). Supprimer ce dossier uniquement si aucun pipeline ne tourne."
fi
cleanup_lock() {
    rmdir "$LOCK_DIR" 2>/dev/null || true
}
trap cleanup_lock EXIT

export NXF_HOME="${JEDI_ROOT}/.nextflow"
export NXF_OPTS="${NXF_OPTS:--Xms1g -Xmx4g}"
export NXF_SINGULARITY_CACHEDIR="$CONTAINER_CACHE"
export NXF_APPTAINER_CACHEDIR="$CONTAINER_CACHE"
export APPTAINER_CACHEDIR="$CONTAINER_CACHE"
export SINGULARITY_CACHEDIR="$CONTAINER_CACHE"
export NXF_SINGULARITY_PULL_TIMEOUT="2h"
export TMPDIR="${TMP_ROOT}/general"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
export XDG_CONFIG_HOME="${TMP_ROOT}/xdg"
export MPLCONFIGDIR="${TMP_ROOT}/mpl"
export NUMBA_CACHE_DIR="${TMP_ROOT}/numba"

command -v nextflow >/dev/null 2>&1 || die "Nextflow est absent du PATH."
[[ -d "$RAW_DIR" ]] || die "Dossier FASTQ absent : ${RAW_DIR}"
[[ -s "$INFO_XLSX" ]] || die "Tableur absent ou vide : ${INFO_XLSX}"

PROFILE="$(choose_profile)" || die "Apptainer, Singularity ou Docker est requis."

for value in "$TRUNCLEN_F" "$TRUNCLEN_R" "$TRUNC_QMIN" "$MAX_EE" "$MIN_LEN" "$RAREFACTION_POINTS"; do
    is_uint "$value" || die "Paramètre entier invalide : ${value}"
done
[[ "$RAREFACTION_POINTS" -ge 2 ]] || die "RAREFACTION_POINTS doit être >= 2."

cat > "${INPUT_DIR}/nextflow_infrastructure.config" <<EOF
singularity {
    autoMounts = true
    cacheDir = '${CONTAINER_CACHE}'
    pullTimeout = '2h'
}
apptainer {
    autoMounts = true
    cacheDir = '${CONTAINER_CACHE}'
    pullTimeout = '2h'
}
process {
    errorStrategy = { task.exitStatus in [137, 140, 143] ? 'retry' : 'terminate' }
    maxRetries = 2
}
EOF

log "Pipeline JEDI corrigé"
log "Racine       : ${JEDI_ROOT}"
log "nf-core      : ampliseq ${NFCORE_VERSION}"
log "Profil       : ${PROFILE}"
log "Classifieurs : ${SILVA_DB} + ${PR2_DB}"
log "Amorces      : ${FW_PRIMER} / ${RV_PRIMER}"
nextflow -version

# -------------------------- CHOIX DU PYTHON EXCEL ----------------------------
PYTHON_CMD=()
if command -v python3 >/dev/null 2>&1 && python3 -c 'import pandas, openpyxl' >/dev/null 2>&1; then
    PYTHON_CMD=(python3)
elif command -v python >/dev/null 2>&1 && python -c 'import pandas, openpyxl' >/dev/null 2>&1; then
    PYTHON_CMD=(python)
elif command -v conda >/dev/null 2>&1 && conda run -n "$EXCEL_ENV" python -c 'import pandas, openpyxl' >/dev/null 2>&1; then
    PYTHON_CMD=(conda run --no-capture-output -n "$EXCEL_ENV" python)
else
    die "Aucun Python avec pandas+openpyxl. Installer ces modules ou vérifier l'environnement Conda ${EXCEL_ENV}."
fi

SAMPLESHEET="${INPUT_DIR}/samplesheet_jedi.tsv"
METADATA="${INPUT_DIR}/metadata_jedi.tsv"
ID_MAP="${INPUT_DIR}/sample_id_mapping.tsv"
CONTROL_FLAG="${INPUT_DIR}/controls_detected.flag"
rm -f "$CONTROL_FLAG"

log "Construction et validation du samplesheet depuis ${INFO_XLSX}"
"${PYTHON_CMD[@]}" - "$INFO_XLSX" "$RAW_DIR" "$SAMPLESHEET" "$METADATA" "$ID_MAP" "$CONTROL_FLAG" <<'PY'
import re
import sys
import unicodedata
from pathlib import Path
import pandas as pd

xlsx, raw_dir, samplesheet, metadata, id_map, control_flag = sys.argv[1:]
raw_dir = Path(raw_dir).resolve()
df = pd.read_excel(xlsx, dtype=str).fillna("")
df.columns = [str(c).strip() for c in df.columns]

required = ["R1", "R2", "New_label"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit("Colonnes obligatoires absentes : " + ", ".join(missing))
for c in df.columns:
    df[c] = df[c].astype(str).str.strip()

used = (df["R1"] != "") | (df["R2"] != "") | (df["New_label"] != "")
partial = df.loc[used & ((df["R1"] == "") | (df["R2"] == "") | (df["New_label"] == "")), required]
if not partial.empty:
    raise SystemExit("Lignes incomplètes R1/R2/New_label :\n" + partial.to_string(index=False))
df = df.loc[(df["R1"] != "") & (df["R2"] != "") & (df["New_label"] != "")].copy()
if df.empty:
    raise SystemExit("Aucun échantillon complet dans le tableur.")

def ascii_text(value):
    return unicodedata.normalize("NFKD", str(value)).encode("ascii", "ignore").decode()

def make_sample_id(value):
    value = re.sub(r"[^A-Za-z0-9_]+", "_", ascii_text(value).strip())
    value = re.sub(r"_+", "_", value).strip("_")
    if not value:
        raise SystemExit("Un New_label ne produit aucun identifiant valide.")
    if not value[0].isalpha():
        value = "S_" + value
    return value[:36]

df["sample"] = [make_sample_id(x) for x in df["New_label"]]
if df["sample"].duplicated().any():
    dup = df.loc[df["sample"].duplicated(keep=False), ["New_label", "sample"]]
    raise SystemExit("Identifiants dupliqués après normalisation :\n" + dup.to_string(index=False))

def fastq_path(value):
    p = Path(value)
    if not p.is_absolute():
        p = raw_dir / p
    p = p.resolve()
    if not p.is_file() or p.stat().st_size == 0:
        raise SystemExit(f"FASTQ absent ou vide : {p}")
    if not (str(p).endswith(".fastq.gz") or str(p).endswith(".fq.gz")):
        raise SystemExit(f"FASTQ non compressé ou extension invalide : {p}")
    return str(p)

out = pd.DataFrame({
    "sample": df["sample"],
    "fastq_1": [fastq_path(x) for x in df["R1"]],
    "fastq_2": [fastq_path(x) for x in df["R2"]],
})

def normalized_name(value):
    return re.sub(r"[^a-z0-9]+", "", ascii_text(value).lower())

by_norm = {normalized_name(c): c for c in df.columns}
run_col = next((by_norm[x] for x in ("run", "sequencingrun", "runid") if x in by_norm), None)
control_col = next((by_norm[x] for x in ("control", "controle", "negativecontrol") if x in by_norm), None)
quant_col = next((by_norm[x] for x in ("quantreading", "dnaquantity", "dnaconcentration") if x in by_norm), None)

if run_col and (df[run_col] != "").any():
    out["run"] = df[run_col].replace("", "run1")

if control_col:
    def parse_control(value):
        x = ascii_text(value).strip().lower()
        controls = {"control", "controle", "negative", "negatif", "blank", "blanc", "ntc"}
        samples = {"", "sample", "echantillon", "positive", "positif"}
        if x in controls or "negative control" in x or "controle negatif" in x:
            return "control"
        if x in samples:
            return "sample"
        raise SystemExit(f"Valeur de contrôle non interprétable : {value!r}")
    out["control"] = [parse_control(x) for x in df[control_col]]
    if (out["control"] == "control").any():
        Path(control_flag).write_text("controls present\n")

if quant_col and (df[quant_col] != "").all():
    q = pd.to_numeric(df[quant_col].str.replace(",", ".", regex=False), errors="coerce")
    if q.notna().all():
        out["quant_reading"] = q

out.to_csv(samplesheet, sep="\t", index=False, lineterminator="\n")
pd.DataFrame({"New_label_original": df["New_label"], "sample": df["sample"]}).to_csv(
    id_map, sep="\t", index=False, lineterminator="\n"
)

meta = df.drop(columns=["R1", "R2", "New_label", "sample"], errors="ignore").copy()
seen = {}
new_columns = []
for col in meta.columns:
    base = re.sub(r"[^A-Za-z0-9_]+", "_", ascii_text(col)).strip("_") or "metadata"
    seen[base] = seen.get(base, 0) + 1
    new_columns.append(base if seen[base] == 1 else f"{base}_{seen[base]}")
meta.columns = new_columns
meta.insert(0, "ID", df["sample"].values)
meta.to_csv(metadata, sep="\t", index=False, lineterminator="\n")
print(f"{len(out)} échantillons écrits dans {samplesheet}")
PY

[[ -s "$SAMPLESHEET" ]] || die "Le samplesheet n'a pas été produit."
[[ -s "$METADATA" ]] || die "Le fichier de métadonnées n'a pas été produit."

log "Contrôle gzip de tous les FASTQ"
while IFS=$'\t' read -r sample fastq_1 fastq_2 rest; do
    [[ "$sample" == "sample" ]] && continue
    gzip -t "$fastq_1"
    gzip -t "$fastq_2"
done < "$SAMPLESHEET"

cp -f "$0" "${INPUT_DIR}/pipeline_jedi_fixed_executed.sh" 2>/dev/null || true

# ----------------------------- PARAMÈTRES SILVA ------------------------------
SILVA_PARAMS="${INPUT_DIR}/params_silva_fixed.yaml"
cat > "$SILVA_PARAMS" <<YAML
input: "${SAMPLESHEET}"
outdir: "${SILVA_OUT}"
FW_primer: "${FW_PRIMER}"
RV_primer: "${RV_PRIMER}"
ref_taxonomy_storage: "${DB_CACHE}"
save_intermediates: true

mergepairs_strategy: "consensus"
mergepairs_consensus_match: 1
mergepairs_consensus_mismatch: -2
mergepairs_consensus_gap: -4
mergepairs_consensus_minoverlap: 12
mergepairs_consensus_maxmismatch: 0
mergepairs_consensus_percentile_cutoff: 0.001

trunclenf: ${TRUNCLEN_F}
trunclenr: ${TRUNCLEN_R}
trunc_qmin: ${TRUNC_QMIN}
trunc_rmin: ${TRUNC_RMIN}
max_ee: ${MAX_EE}
min_len: ${MIN_LEN}
sample_inference: "independent"

dada_ref_taxonomy: "${SILVA_DB}"
cut_dada_ref_taxonomy: true
skip_dada_addspecies: true
exclude_taxa: "none"
min_frequency: 1
min_samples: 1

# Correctif essentiel : aucune étape QIIME2/Emperor n'est lancée.
skip_qiime: true
skip_qiime_downstream: true
skip_alpha_rarefaction: true
skip_diversity_indices: true
skip_abundance_tables: true
report_title: "JEDI valormicro - ASV et SILVA"
YAML

# --------------------------- INFÉRENCE + SILVA -------------------------------
if [[ "$FORCE_ALL" == "1" || "$FORCE_SILVA" == "1" ]] || ! silva_core_complete; then
    log "ÉTAPE 1/3 : DADA2/JEDI et classification SILVA"
    run_nextflow "$SILVA_PARAMS" "${WORK_DIR}/silva" SILVA
else
    log "ÉTAPE 1/3 : résultats SILVA centraux déjà complets ; réutilisation."
fi

ASV_FASTA="${SILVA_OUT}/dada2/ASV_seqs.fasta"
ASV_TABLE="${SILVA_OUT}/dada2/ASV_table.tsv"
[[ -s "$ASV_FASTA" ]] || die "ASV FASTA absent : ${ASV_FASTA}"
[[ -s "$ASV_TABLE" ]] || die "Table ASV absente : ${ASV_TABLE}"
has_taxonomy "$SILVA_OUT" || die "Taxonomie SILVA absente dans ${SILVA_OUT}/dada2"

# ------------------------------- TAXONOMIE PR2 -------------------------------
PR2_PARAMS="${INPUT_DIR}/params_pr2_fixed.yaml"
cat > "$PR2_PARAMS" <<YAML
input_fasta: "${ASV_FASTA}"
outdir: "${PR2_OUT}"
FW_primer: "${FW_PRIMER}"
RV_primer: "${RV_PRIMER}"
ref_taxonomy_storage: "${DB_CACHE}"
save_intermediates: true

dada_ref_taxonomy: "${PR2_DB}"
cut_dada_ref_taxonomy: true
skip_dada_addspecies: true

skip_fastqc: true
skip_barrnap: true
skip_qiime: true
skip_qiime_downstream: true
skip_alpha_rarefaction: true
skip_diversity_indices: true
skip_abundance_tables: true
report_title: "JEDI valormicro - classification PR2"
YAML

if [[ "$FORCE_ALL" == "1" || "$FORCE_PR2" == "1" ]] || ! pr2_core_complete; then
    log "ÉTAPE 2/3 : classification des mêmes ASV avec PR2"
    run_nextflow "$PR2_PARAMS" "${WORK_DIR}/pr2" PR2
else
    log "ÉTAPE 2/3 : taxonomie PR2 déjà complète ; réutilisation."
fi
has_taxonomy "$PR2_OUT" || die "Taxonomie PR2 absente dans ${PR2_OUT}/dada2"

# -------------------------- INTÉGRATION CROSS-DOMAIN -------------------------
log "ÉTAPE 3/3 : intégration SILVA/PR2 et analyses tabulaires"
"${PYTHON_CMD[@]}" - "$SILVA_OUT" "$PR2_OUT" "$ASV_TABLE" "$ASV_FASTA" \
    "$INTEGRATED_OUT" "$RAREFACTION_POINTS" <<'PY'
import csv
import math
import shutil
import sys
from collections import defaultdict
from pathlib import Path

silva_root, pr2_root, table_path, fasta_path, outdir, npoints = sys.argv[1:]
silva_root = Path(silva_root)
pr2_root = Path(pr2_root)
table_path = Path(table_path)
fasta_path = Path(fasta_path)
outdir = Path(outdir)
npoints = int(npoints)
outdir.mkdir(parents=True, exist_ok=True)
(outdir / "tables_by_domain").mkdir(exist_ok=True)

def taxonomy_file(root, wanted):
    candidates = []
    for p in (root / "dada2").glob("ASV_tax.*.tsv"):
        low = p.name.lower()
        if p.is_file() and p.stat().st_size and "species" not in low and "into-qiime" not in low:
            candidates.append(p)
    if not candidates:
        raise SystemExit(f"Aucune taxonomie dans {root / 'dada2'}")
    candidates.sort(key=lambda p: (wanted.lower() not in p.name.lower(), len(p.name), p.name))
    return candidates[0]

def read_taxonomy(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"Taxonomie sans en-tête : {path}")
        id_candidates = {"asv_id", "feature id", "featureid", "id", "#otu id"}
        id_col = next((c for c in reader.fieldnames if c.strip().lower() in id_candidates), reader.fieldnames[0])
        ignored = {"sequence", "confidence", "database"}
        rank_cols = [c for c in reader.fieldnames if c != id_col and c.strip().lower() not in ignored]
        data = {}
        for row in reader:
            asv = (row.get(id_col) or "").strip()
            if asv:
                data[asv] = [(row.get(c) or "").strip() for c in rank_cols]
    return data

def read_fasta(path):
    seqs, current = {}, None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in seqs:
                    raise SystemExit(f"Identifiant FASTA dupliqué : {current}")
                seqs[current] = []
            elif current is None:
                raise SystemExit("FASTA invalide : séquence avant le premier en-tête")
            else:
                seqs[current].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}

def clean_taxon(value):
    value = str(value).strip().strip(";")
    return value

def informative(ranks):
    bad = {"", "na", "nan", "none", "unassigned", "unknown", "root"}
    return [x for x in (clean_taxon(v) for v in ranks) if x.lower() not in bad]

def tax_string(ranks):
    return ";".join(informative(ranks)) or "Unassigned"

def contains(text, terms):
    low = text.lower()
    return any(term in low for term in terms)

def safe_name(value):
    return "".join(c if c.isalnum() or c in "_-" else "_" for c in value)

silva_path = taxonomy_file(silva_root, "silva")
pr2_path = taxonomy_file(pr2_root, "pr2")
silva = read_taxonomy(silva_path)
pr2 = read_taxonomy(pr2_path)
seqs = read_fasta(fasta_path)

with table_path.open(encoding="utf-8-sig", newline="") as handle:
    raw_rows = list(csv.reader(handle, delimiter="\t"))
if len(raw_rows) < 2:
    raise SystemExit(f"Table ASV vide ou invalide : {table_path}")

header = raw_rows[0]
if len(header) < 2:
    raise SystemExit("La table ASV ne contient aucune colonne échantillon.")
sample_names = header[1:]
if len(set(sample_names)) != len(sample_names):
    raise SystemExit("Noms d'échantillons dupliqués dans la table ASV.")

counts = {}
for line_no, row in enumerate(raw_rows[1:], 2):
    if not row or not row[0].strip():
        continue
    if len(row) != len(header):
        raise SystemExit(f"Nombre de colonnes incorrect ligne {line_no} de {table_path}")
    asv = row[0].strip()
    if asv in counts:
        raise SystemExit(f"ASV dupliqué dans la table : {asv}")
    try:
        vals = [int(float(x or 0)) for x in row[1:]]
    except ValueError as exc:
        raise SystemExit(f"Comptage non numérique ligne {line_no}: {exc}")
    if any(x < 0 for x in vals):
        raise SystemExit(f"Comptage négatif ligne {line_no}")
    counts[asv] = vals

missing_fasta = sorted(set(counts) - set(seqs))
if missing_fasta:
    raise SystemExit(f"{len(missing_fasta)} ASV de la table sont absents du FASTA")

records = {}
for asv in counts:
    s_tax = tax_string(silva.get(asv, []))
    p_tax = tax_string(pr2.get(asv, []))
    sl, pl = s_tax.lower(), p_tax.lower()

    if contains(sl, ["chloroplast", "plastid"]):
        domain, source, chosen = "Eukaryota_plastid", "SILVA", s_tax
    elif contains(sl, ["mitochond"]):
        domain, source, chosen = "Eukaryota_mitochondria", "SILVA", s_tax
    elif contains(sl, ["archaea", "d_0__archaea", "k__archaea"]):
        domain, source, chosen = "Archaea", "SILVA", s_tax
    elif contains(sl, ["bacteria", "d_0__bacteria", "k__bacteria"]):
        domain, source, chosen = "Bacteria", "SILVA", s_tax
    elif contains(pl, ["eukaryota", "eukaryote", "kingdom__eukaryota"]):
        domain, source, chosen = "Eukaryota", "PR2", p_tax
    elif p_tax != "Unassigned":
        domain, source, chosen = "Eukaryota_candidate", "PR2", p_tax
    elif s_tax != "Unassigned":
        domain, source, chosen = "Unresolved_SILVA", "SILVA", s_tax
    else:
        domain, source, chosen = "Unassigned", "none", "Unassigned"

    records[asv] = {
        "domain": domain,
        "source": source,
        "chosen": chosen,
        "silva": s_tax,
        "pr2": p_tax,
        "length": len(seqs[asv]),
    }

ordered_asvs = list(counts)
domains = sorted({records[a]["domain"] for a in ordered_asvs})

# Taxonomie intégrée.
with (outdir / "taxonomy_JEDI_consensus.tsv").open("w", encoding="utf-8", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["ASV_ID", "JEDI_domain", "selected_source", "selected_taxonomy",
                "SILVA_taxonomy", "PR2_taxonomy", "ASV_length"])
    for asv in ordered_asvs:
        r = records[asv]
        w.writerow([asv, r["domain"], r["source"], r["chosen"], r["silva"], r["pr2"], r["length"]])

# Tables ASV brute et enrichie.
with (outdir / "ASV_table_JEDI_counts.tsv").open("w", encoding="utf-8", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["ASV_ID"] + sample_names)
    for asv in ordered_asvs:
        w.writerow([asv] + counts[asv])

with (outdir / "ASV_table_JEDI_with_taxonomy.tsv").open("w", encoding="utf-8", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["ASV_ID", "JEDI_domain", "selected_source", "selected_taxonomy",
                "SILVA_taxonomy", "PR2_taxonomy", "ASV_length"] + sample_names)
    for asv in ordered_asvs:
        r = records[asv]
        w.writerow([asv, r["domain"], r["source"], r["chosen"], r["silva"], r["pr2"], r["length"]] + counts[asv])

# Une table de comptages par domaine, sans supprimer les tables globales.
for domain in domains:
    path = outdir / "tables_by_domain" / f"ASV_counts_{safe_name(domain)}.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["ASV_ID", "selected_taxonomy"] + sample_names)
        for asv in ordered_asvs:
            if records[asv]["domain"] == domain:
                w.writerow([asv, records[asv]["chosen"]] + counts[asv])

# Synthèses par domaine.
domain_counts = {d: [0] * len(sample_names) for d in domains}
domain_richness = defaultdict(int)
for asv in ordered_asvs:
    d = records[asv]["domain"]
    domain_richness[d] += 1
    domain_counts[d] = [a + b for a, b in zip(domain_counts[d], counts[asv])]

with (outdir / "domain_counts_per_sample.tsv").open("w", encoding="utf-8", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["JEDI_domain"] + sample_names)
    for d in domains:
        w.writerow([d] + domain_counts[d])

sample_totals = [sum(counts[a][i] for a in ordered_asvs) for i in range(len(sample_names))]
with (outdir / "domain_relative_abundance_per_sample.tsv").open("w", encoding="utf-8", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["JEDI_domain"] + sample_names)
    for d in domains:
        rel = [domain_counts[d][i] / sample_totals[i] if sample_totals[i] else 0.0 for i in range(len(sample_names))]
        w.writerow([d] + [f"{x:.10f}" for x in rel])

with (outdir / "domain_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["JEDI_domain", "ASV_richness", "total_reads"])
    for d in domains:
        w.writerow([d, domain_richness[d], sum(domain_counts[d])])

# Alpha-diversité non raréfiée.
with (outdir / "alpha_diversity_unrarefied.tsv").open("w", encoding="utf-8", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["sample", "total_reads", "observed_ASVs", "Shannon", "Simpson_1_minus_D"])
    for i, sample in enumerate(sample_names):
        vals = [counts[a][i] for a in ordered_asvs if counts[a][i] > 0]
        total = sum(vals)
        observed = len(vals)
        if total:
            proportions = [x / total for x in vals]
            shannon = -sum(p * math.log(p) for p in proportions)
            simpson = 1.0 - sum(p * p for p in proportions)
        else:
            shannon = simpson = 0.0
        w.writerow([sample, total, observed, f"{shannon:.10f}", f"{simpson:.10f}"])

# Raréfaction analytique : espérance du nombre d'ASV observés à profondeur n.
def logchoose(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def expected_richness(abundances, depth):
    total = sum(abundances)
    if depth <= 0 or total == 0:
        return 0.0
    if depth >= total:
        return float(sum(x > 0 for x in abundances))
    denominator = logchoose(total, depth)
    expected = 0.0
    for abundance in abundances:
        if abundance <= 0:
            continue
        if total - abundance < depth:
            p_absent = 0.0
        else:
            p_absent = math.exp(logchoose(total - abundance, depth) - denominator)
        expected += 1.0 - p_absent
    return expected

with (outdir / "rarefaction_expected_ASVs.tsv").open("w", encoding="utf-8", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["sample", "depth", "expected_observed_ASVs"])
    for i, sample in enumerate(sample_names):
        abund = [counts[a][i] for a in ordered_asvs]
        total = sum(abund)
        if total == 0:
            w.writerow([sample, 0, "0.000000"])
            continue
        depths = {1, total}
        for j in range(1, npoints + 1):
            depths.add(max(1, round(total * j / npoints)))
        for depth in sorted(depths):
            w.writerow([sample, depth, f"{expected_richness(abund, depth):.6f}"])

shutil.copy2(fasta_path, outdir / "ASV_sequences_JEDI.fasta")
shutil.copy2(silva_path, outdir / "taxonomy_SILVA_original.tsv")
shutil.copy2(pr2_path, outdir / "taxonomy_PR2_original.tsv")

print(f"Taxonomie SILVA : {silva_path}")
print(f"Taxonomie PR2   : {pr2_path}")
print(f"ASV intégrés    : {len(ordered_asvs)}")
print(f"Échantillons    : {len(sample_names)}")
print("Domaines        : " + ", ".join(domains))
PY

# ---------------------------- CONTRÔLES FINAUX -------------------------------
REQUIRED_OUTPUTS=(
    "${INTEGRATED_OUT}/ASV_table_JEDI_counts.tsv"
    "${INTEGRATED_OUT}/ASV_table_JEDI_with_taxonomy.tsv"
    "${INTEGRATED_OUT}/taxonomy_JEDI_consensus.tsv"
    "${INTEGRATED_OUT}/domain_summary.tsv"
    "${INTEGRATED_OUT}/domain_counts_per_sample.tsv"
    "${INTEGRATED_OUT}/domain_relative_abundance_per_sample.tsv"
    "${INTEGRATED_OUT}/alpha_diversity_unrarefied.tsv"
    "${INTEGRATED_OUT}/rarefaction_expected_ASVs.tsv"
    "${INTEGRATED_OUT}/ASV_sequences_JEDI.fasta"
)
for output in "${REQUIRED_OUTPUTS[@]}"; do
    [[ -s "$output" ]] || die "Sortie finale absente ou vide : ${output}"
done

{
    printf 'parameter\tvalue\n'
    printf 'date\t%s\n' "$(date --iso-8601=seconds)"
    printf 'project_dir\t%s\n' "$PROJECT_DIR"
    printf 'jedi_root\t%s\n' "$JEDI_ROOT"
    printf 'nfcore_ampliseq\t%s\n' "$NFCORE_VERSION"
    printf 'profile\t%s\n' "$PROFILE"
    printf 'forward_primer\t%s\n' "$FW_PRIMER"
    printf 'reverse_primer\t%s\n' "$RV_PRIMER"
    printf 'silva_database\t%s\n' "$SILVA_DB"
    printf 'pr2_database\t%s\n' "$PR2_DB"
    printf 'mergepairs_strategy\tconsensus\n'
    printf 'trunclen_f\t%s\n' "$TRUNCLEN_F"
    printf 'trunclen_r\t%s\n' "$TRUNCLEN_R"
    printf 'max_ee\t%s\n' "$MAX_EE"
    printf 'qiime_downstream\tdisabled_to_avoid_emperor_tmp_failure\n'
} > "${INTEGRATED_OUT}/run_manifest.tsv"

# Marqueur écrit uniquement après validation de toutes les sorties.
printf 'SUCCESS\t%s\n' "$(date --iso-8601=seconds)" > "${INTEGRATED_OUT}/PIPELINE_SUCCESS.txt"

log "Pipeline JEDI terminé avec succès."
log "Résultats intégrés : ${INTEGRATED_OUT}"
log "Table principale   : ${INTEGRATED_OUT}/ASV_table_JEDI_with_taxonomy.tsv"
log "Résumé domaines    : ${INTEGRATED_OUT}/domain_summary.tsv"
log "Raréfaction        : ${INTEGRATED_OUT}/rarefaction_expected_ASVs.tsv"
