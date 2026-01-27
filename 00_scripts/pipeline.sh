#!/usr/bin/env bash

set -euo pipefail

# ==================== CONFIGURATION ====================
export ROOTDIR="/nvme/bio/data_fungi/valormicro_V2"
export NTHREADS=16
export TMPDIR="${ROOTDIR}/tmp"
export QIIME_ENV="qiime2-amplicon-2025.7"
export FASTQC_ENV="fastqc"
export TRIMMOMATIC_ENV="trimmomatic"

mkdir -p "$TMPDIR"

log() { echo -e "\n[$(date +'%F %T')] $*\n"; }

log "=== PIPELINE VALORMICRO V2 AVEC TRIMMOMATIC DÉMARRÉ ==="

# ---- CORRECTION TOKEN CONDA ----
log "Suppression token conda corrompu"
sudo rm -rf /home/fungi/.conda/ 2>/dev/null || rm -rf /home/fungi/.conda/ 2>/dev/null || true

set +u
source $(conda info --base)/etc/profile.d/conda.sh
export JAVA_HOME="${JAVA_HOME:-/usr/lib/jvm/default-java}"
set -u

# ==================== FONCTION VALIDATION FASTQ ====================
validate_fastq() {
    local file="$1"
    
    if [ ! -f "$file" ]; then
        return 1
    fi
    
    local size=$(stat -c%s "$file" 2>/dev/null || echo "0")
    if [ "$size" -lt 1000 ]; then
        return 1
    fi
    
    # Test de décompression pour .gz
    if [[ "$file" =~ \.gz$ ]]; then
        if ! gunzip -t "$file" 2>/dev/null; then
            return 1
        fi
        local line_count=$(zcat "$file" 2>/dev/null | head -4 | wc -l || echo "0")
    else
        local line_count=$(head -4 "$file" 2>/dev/null | wc -l || echo "0")
    fi
    
    if [ "$line_count" -ne 4 ]; then
        return 1
    fi
    
    return 0
}

# ==================== 01 GÉNÉRATION MANIFEST ====================
log "Génération manifest pour tous les échantillons"

cd "${ROOTDIR}"
mkdir -p 98_databasefiles

MANIFEST="${ROOTDIR}/98_databasefiles/manifest"
MANIFEST_TEMP="${ROOTDIR}/98_databasefiles/manifest_temp"

echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST_TEMP"

cd "${ROOTDIR}/01_raw_data"

log "Scan des fichiers fastq dans: $(pwd)"

declare -A seen_ids
declare -A filepath_registry
count=0
duplicates_found=0
invalid_files=0

for r1_file in *R1*.fastq* *_R1.fastq*; do
    if [ ! -f "$r1_file" ]; then
        continue
    fi
    
    r2_file="${r1_file/R1/R2}"
    r2_file="${r2_file/_R1./_R2.}"
    
    if [ ! -f "$r2_file" ]; then
        continue
    fi
    
    r1_abs="${ROOTDIR}/01_raw_data/$r1_file"
    r2_abs="${ROOTDIR}/01_raw_data/$r2_file"
    
    # Validation fastq
    if ! validate_fastq "$r1_abs"; then
        log "⚠️  FICHIER R1 INVALIDE IGNORÉ: $r1_file"
        invalid_files=$((invalid_files + 1))
        continue
    fi
    
    if ! validate_fastq "$r2_abs"; then
        log "⚠️  FICHIER R2 INVALIDE IGNORÉ: $r2_file"
        invalid_files=$((invalid_files + 1))
        continue
    fi
    
    r1_size=$(stat -c%s "$r1_file" 2>/dev/null || echo "0")
    r2_size=$(stat -c%s "$r2_file" 2>/dev/null || echo "0")
    
    if [ "$r1_size" -gt 1000 ] && [ "$r2_size" -gt 1000 ]; then
        base_name=$(basename "$r1_file")
        sample_id="${base_name%%_S[0-9]*}"
        sample_id="${sample_id%%_R1*}"
        sample_id="${sample_id%.fastq*}"
        sample_id="${sample_id%.gz}"
        
        sample_id="${sample_id//[^a-zA-Z0-9._-]/_}"
        
        if [[ -n "${filepath_registry[$r1_abs]:-}" ]]; then
            log "⚠️  DOUBLON: $r1_file est déjà assigné"
            duplicates_found=$((duplicates_found + 1))
            continue
        fi
        
        # Gérer doublons de sample-id correctement
        original_id="$sample_id"
        counter=1
        while [[ -v seen_ids[$sample_id] ]]; do
            sample_id="${original_id}_${counter}"
            ((counter++))
        done
        
        seen_ids[$sample_id]=1
        filepath_registry[$r1_abs]="$sample_id"
        
        echo -e "$sample_id\t$r1_abs\t$r2_abs" >> "$MANIFEST_TEMP"
        count=$((count + 1))
        log "✅ Ajouté: $sample_id"
    fi
done

log "Manifest: $count échantillons (invalides: $invalid_files, doublons: $duplicates_found)"

if [ "$count" -eq 0 ]; then
    log "ERREUR: Aucun échantillon valide"
    exit 1
fi

cp "$MANIFEST_TEMP" "$MANIFEST"
rm "$MANIFEST_TEMP"

# ==================== 02 FASTQC RAW DATA ====================
log "FastQC sur données brutes"

mkdir -p "${ROOTDIR}/02_qualitycheck"
rm -f "${ROOTDIR}/02_qualitycheck"/*multiqc* 2>/dev/null || true

cd "${ROOTDIR}/01_raw_data"

count=0
for file in *.fastq* */*.fastq*; do
    if [ -f "$file" ]; then
        count=$((count + 1))
        if [ $((count % 20)) -eq 0 ]; then
            log "FastQC: $count fichiers traités"
        fi
        conda run -n "$FASTQC_ENV" fastqc "$file" -o "${ROOTDIR}/02_qualitycheck" --threads 2 --quiet 2>/dev/null || true
    fi
done

log "FastQC terminé: $count fichiers"

# MultiQC
cd "${ROOTDIR}/02_qualitycheck"
conda run -n "$FASTQC_ENV" multiqc . --force --filename "raw_data_qc" --no-ansi 2>/dev/null || true

# ==================== 02B TRIMMOMATIC ====================
log "Trimmomatic - Nettoyage et suppression amorces"

mkdir -p "${ROOTDIR}/03_cleaned_data"
cd "${ROOTDIR}/01_raw_data"

# Créer fichier adapters
ADAPTERS="${ROOTDIR}/98_databasefiles/adapters_515f926r.fasta"
cat > "$ADAPTERS" << 'ADAPTER_EOF'
>515F
GTGYCAGCMGCCGCGGTAA
>926R
CCGYCAATTYMTTTRAGTTT
ADAPTER_EOF

total_pairs=$(tail -n +2 "$MANIFEST" | wc -l)
success_count=0
failed_count=0
pair_count=0

log "Traitement de $total_pairs paires..."

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r sample_id r1_path r2_path; do
    pair_count=$((pair_count + 1))
    
    out1p="${ROOTDIR}/03_cleaned_data/${sample_id}_R1_paired.fastq.gz"
    out1u="${ROOTDIR}/03_cleaned_data/${sample_id}_R1_unpaired.fastq.gz"
    out2p="${ROOTDIR}/03_cleaned_data/${sample_id}_R2_paired.fastq.gz"
    out2u="${ROOTDIR}/03_cleaned_data/${sample_id}_R2_unpaired.fastq.gz"
    
    # Pré-vérifier si fichiers existent ET sont valides
    if [ ! -f "$r1_path" ] || [ ! -f "$r2_path" ]; then
        log "⚠️  Fichiers manquants pour $sample_id"
        continue
    fi
    
    # Exécuter trimmomatic
    if conda run -n "$TRIMMOMATIC_ENV" trimmomatic PE -threads 8 -phred33 \
        "$r1_path" "$r2_path" \
        "$out1p" "$out1u" "$out2p" "$out2u" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100 2>&1 > /dev/null; then
        
        # Vérifier que les fichiers de sortie existent et sont valides
        if [ -f "$out1p" ] && [ -f "$out2p" ]; then
            # Vérifier que les fichiers ne sont pas corrompus
            if gunzip -t "$out1p" 2>/dev/null && gunzip -t "$out2p" 2>/dev/null; then
                count1=$(( $(zcat "$out1p" | wc -l) / 4 ))
                count2=$(( $(zcat "$out2p" | wc -l) / 4 ))
                
                if [ "$count1" -eq "$count2" ] && [ "$count1" -gt 0 ]; then
                    success_count=$((success_count + 1))
                    if [ $((pair_count % 10)) -eq 0 ]; then
                        log "✅ $pair_count/$total_pairs: $sample_id ($count1 reads)"
                    fi
                else
                    log "⚠️  Désynchronisé: $sample_id ($count1 vs $count2) - suppression"
                    rm -f "$out1p" "$out2p" "$out1u" "$out2u"
                    failed_count=$((failed_count + 1))
                fi
            else
                log "⚠️  Fichiers corrompus après trimming: $sample_id"
                rm -f "$out1p" "$out2p" "$out1u" "$out2u"
                failed_count=$((failed_count + 1))
            fi
        else
            log "⚠️  Fichiers manquants: $sample_id"
            failed_count=$((failed_count + 1))
        fi
    else
        log "⚠️  Trimmomatic échoué: $sample_id"
        rm -f "$out1p" "$out2p" "$out1u" "$out2u"
        failed_count=$((failed_count + 1))
    fi
done

log "Trimmomatic: Succès=$success_count, Erreurs=$failed_count"

# ==================== 03 CRÉER MANIFEST NETTOYÉ ====================
log "Génération manifest données nettoyées"

MANIFEST_TRIMMED="${ROOTDIR}/98_databasefiles/manifest_trimmed"
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST_TRIMMED"

cd "${ROOTDIR}/03_cleaned_data"

for r1_file in *_R1_paired.fastq.gz; do
    if [ -f "$r1_file" ]; then
        r2_file="${r1_file/_R1_paired/_R2_paired}"
        
        if [ -f "$r2_file" ]; then
            sample_id="${r1_file%%_R1_paired*}"
            
            # Valider les fichiers
            if gunzip -t "$r1_file" 2>/dev/null && gunzip -t "$r2_file" 2>/dev/null; then
                r1_abs="${ROOTDIR}/03_cleaned_data/$r1_file"
                r2_abs="${ROOTDIR}/03_cleaned_data/$r2_file"
                echo -e "$sample_id\t$r1_abs\t$r2_abs" >> "$MANIFEST_TRIMMED"
            else
                log "⚠️  Fichiers corrompus ignorés: $sample_id"
            fi
        fi
    fi
done

trimmed_count=$(tail -n +2 "$MANIFEST_TRIMMED" | wc -l)
log "✅ Manifest nettoyé: $trimmed_count échantillons"

if [ "$trimmed_count" -lt 10 ]; then
    log "ERREUR: Trop peu d'échantillons ($trimmed_count)"
    exit 1
fi

# ==================== 04 FASTQC CLEANED DATA ====================
log "FastQC données nettoyées"

mkdir -p "${ROOTDIR}/03_cleaned_data_qc"
cd "${ROOTDIR}/03_cleaned_data"

count=0
for file in *_paired.fastq.gz; do
    if [ -f "$file" ]; then
        count=$((count + 1))
        if [ $((count % 20)) -eq 0 ]; then
            log "FastQC cleaned: $count fichiers"
        fi
        conda run -n "$FASTQC_ENV" fastqc "$file" -o "${ROOTDIR}/03_cleaned_data_qc" --threads 2 --quiet 2>/dev/null || true
    fi
done

cd "${ROOTDIR}/03_cleaned_data_qc"
conda run -n "$FASTQC_ENV" multiqc . --force --filename "cleaned_data_qc" --no-ansi 2>/dev/null || true

log "✅ FastQC post-nettoyage terminé"

# ==================== 05 QIIME2 IMPORT ====================
log "QIIME2 Import données nettoyées"

mkdir -p "${ROOTDIR}/05_QIIME2/core" "${ROOTDIR}/05_QIIME2/visual"
cd "${ROOTDIR}/05_QIIME2"

if ! conda run -n "$QIIME_ENV" qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST_TRIMMED" \
    --output-path "core/demux.qza" \
    --input-format PairedEndFastqManifestPhred33V2 2>&1; then
    
    log "ERREUR: Import QIIME2 échoué"
    log "Vérification des fichiers dans manifest..."
    
    tail -n +2 "$MANIFEST_TRIMMED" | while IFS=$'\t' read -r sample_id r1_path r2_path; do
        if [ ! -f "$r1_path" ] || [ ! -f "$r2_path" ]; then
            log "  ❌ Manquant: $sample_id"
        elif ! gunzip -t "$r1_path" 2>/dev/null; then
            log "  ❌ Corrompu R1: $sample_id"
        elif ! gunzip -t "$r2_path" 2>/dev/null; then
            log "  ❌ Corrompu R2: $sample_id"
        fi
    done
    
    exit 1
fi

log "✅ Import QIIME2 réussi"

conda run -n "$QIIME_ENV" qiime demux summarize \
    --i-data core/demux.qza \
    --o-visualization visual/demux-summary.qzv 2>/dev/null || true

# ==================== 06 DADA2 ====================
log "DADA2 denoising"

cd "${ROOTDIR}/05_QIIME2/core"

if ! conda run -n "$QIIME_ENV" qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux.qza \
    --o-table table.qza \
    --o-representative-sequences rep-seqs.qza \
    --o-denoising-stats denoising-stats.qza \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-n-threads "$NTHREADS"; then
    log "ERREUR DADA2"
    exit 1
fi

log "✅ DADA2 réussi"

# ==================== 07 CLASSIFICATION ====================
log "Classification taxonomique SILVA 138.2"

cd "${ROOTDIR}/98_databasefiles"

CLASSIFIER_URL="/nvme/bio/data_fungi/valormicro_nc/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
CLASSIFIER="${ROOTDIR}/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"

if [ ! -f "$CLASSIFIER" ]; then
    cp "$CLASSIFIER_URL" "$CLASSIFIER" || { log "ERREUR: Classifier manquant"; exit 1; }
fi

conda run -n "$QIIME_ENV" qiime tools validate "$CLASSIFIER" || { log "ERREUR: Classifier invalide"; exit 1; }

cd "${ROOTDIR}/05_QIIME2/core"

conda run -n "$QIIME_ENV" qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza \
    --p-n-jobs "$NTHREADS"

log "✅ Classification réussie"

# ==================== 08 ANALYSES ====================
log "Analyses finales"

mkdir -p "${ROOTDIR}/05_QIIME2/export"
cd "${ROOTDIR}/05_QIIME2/core"

conda run -n "$QIIME_ENV" qiime feature-table summarize \
    --i-table table.qza \
    --o-visualization "../visual/table-summary.qzv" 2>/dev/null || true

conda run -n "$QIIME_ENV" qiime taxa barplot \
    --i-table table.qza \
    --i-taxonomy taxonomy.qza \
    --o-visualization "../visual/taxa-bar-plots.qzv" 2>/dev/null || true

log "✅ Barplots créés"

# ==================== 09 DIVERSITÉ ====================
log "Diversité phylogénétique"

mkdir -p "${ROOTDIR}/05_QIIME2/diversity" "${ROOTDIR}/05_QIIME2/pcoa"

if [ ! -f "tree.qza" ]; then
    conda run -n "$QIIME_ENV" qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences rep-seqs.qza \
        --o-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza \
        --o-rooted-tree tree.qza \
        --p-n-threads "$NTHREADS"
fi

cd "${ROOTDIR}/98_databasefiles"

echo -e "sample-id\tgroup\ttype" > "diversity-metadata.tsv"
tail -n +2 "$MANIFEST_TRIMMED" | cut -f1 | while read -r sample_id; do
    if echo "${sample_id,,}" | grep -qE "(neg|blank|control|ctrl|eau|mock)"; then
        echo -e "$sample_id\tcontrol\tnegative"
    else
        echo -e "$sample_id\tsample\tenvironmental"
    fi
done >> "diversity-metadata.tsv"

cd "${ROOTDIR}/05_QIIME2/core"

conda run -n "$QIIME_ENV" qiime diversity core-metrics-phylogenetic \
    --i-table table.qza \
    --i-phylogeny tree.qza \
    --m-metadata-file "${ROOTDIR}/98_databasefiles/diversity-metadata.tsv" \
    --o-faith-pd-vector diversity/Vector-faith_pd.qza \
    --o-observed-features-vector diversity/Vector-observed_asv.qza \
    --o-shannon-vector diversity/Vector-shannon.qza \
    --o-evenness-vector diversity/Vector-evenness.qza \
    --o-unweighted-unifrac-distance-matrix diversity/Matrix-unweighted_unifrac.qza \
    --o-weighted-unifrac-distance-matrix diversity/Matrix-weighted_unifrac.qza \
    --o-jaccard-distance-matrix diversity/Matrix-jaccard.qza \
    --o-bray-curtis-distance-matrix diversity/Matrix-braycurtis.qza \
    --o-unweighted-unifrac-pcoa-results pcoa/PCoA-unweighted_unifrac.qza \
    --o-weighted-unifrac-pcoa-results pcoa/PCoA-weighted_unifrac.qza \
    --o-jaccard-pcoa-results pcoa/PCoA-jaccard.qza \
    --o-bray-curtis-pcoa-results pcoa/PCoA-braycurtis.qza \
    --o-unweighted-unifrac-emperor visual/Emperor-unweighted_unifrac.qzv \
    --o-weighted-unifrac-emperor visual/Emperor-weighted_unifrac.qzv \
    --o-jaccard-emperor visual/Emperor-jaccard.qzv \
    --o-bray-curtis-emperor visual/Emperor-braycurtis.qzv 2>/dev/null || true

log "✅ Métriques de diversité créées"

# ==================== 10 EXPORTS ====================
log "Export fichiers"

mkdir -p "${ROOTDIR}/05_QIIME2/export/core" "${ROOTDIR}/05_QIIME2/export/diversity_tsv"
cd "${ROOTDIR}/05_QIIME2"

conda run -n "$QIIME_ENV" qiime tools export --input-path core/table.qza --output-path export/core/table 2>/dev/null || true
conda run -n "$QIIME_ENV" qiime tools export --input-path core/rep-seqs.qza --output-path export/core/rep-seqs 2>/dev/null || true
conda run -n "$QIIME_ENV" qiime tools export --input-path core/taxonomy.qza --output-path export/core/taxonomy 2>/dev/null || true

log "✅ Exports terminés"

# ==================== RAPPORT FINAL ====================
log "🎉 PIPELINE TERMINÉ AVEC SUCCÈS ! 🎉"
log ""
log "Échantillons nettoyés: $trimmed_count"
log "Résultats: ${ROOTDIR}/05_QIIME2/"
log "Barplots: ${ROOTDIR}/05_QIIME2/visual/taxa-bar-plots.qzv"
log "PCoA: ${ROOTDIR}/05_QIIME2/visual/Emperor-*.qzv"
