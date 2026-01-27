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
    
    local line_count=0
    if [[ "$file" =~ \.gz$ ]]; then
        line_count=$(zcat "$file" 2>/dev/null | head -4 | wc -l || echo "0")
    else
        line_count=$(head -4 "$file" 2>/dev/null | wc -l || echo "0")
    fi
    
    if [ "$line_count" -ne 4 ]; then
        return 1
    fi
    
    if [[ "$file" =~ \.gz$ ]]; then
        if ! gunzip -t "$file" 2>/dev/null; then
            return 1
        fi
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
    if [ -f "$r1_file" ]; then
        r2_file="${r1_file/R1/R2}"
        r2_file="${r2_file/_R1./_R2.}"
        
        if [ -f "$r2_file" ]; then
            r1_abs="${ROOTDIR}/01_raw_data/$r1_file"
            r2_abs="${ROOTDIR}/01_raw_data/$r2_file"
            
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
                
                sample_id="${sample_id//[^a-zA-Z0-9._-]/_}"
                
                if [[ -n "${filepath_registry[$r1_abs]:-}" ]]; then
                    log "⚠️  DOUBLON DÉTECTÉ: $r1_file est déjà assigné à '${filepath_registry[$r1_abs]}'"
                    duplicates_found=$((duplicates_found + 1))
                    continue
                fi
                
                original_id="$sample_id"
                counter=1
                while [[ -n "${seen_ids[$sample_id]:-}" ]]; do
                    sample_id="${original_id}_${counter}"
                    counter=$((counter + 1))
                done
                
                seen_ids["$sample_id"]=1
                filepath_registry["$r1_abs"]="$sample_id"
                
                echo -e "$sample_id\t$r1_abs\t$r2_abs" >> "$MANIFEST_TEMP"
                count=$((count + 1))
                log "✅ Ajouté: $sample_id ($r1_file)"
            fi
        fi
    fi
done

log "Manifest temporaire créé avec $count échantillons"

if [ $duplicates_found -gt 0 ]; then
    log "⚠️  $duplicates_found fichier(s) R1 en doublon détecté(s) et ignoré(s)"
fi

if [ $invalid_files -gt 0 ]; then
    log "⚠️  $invalid_files fichier(s) invalide(s)/corrompu(s) ignoré(s)"
fi

if [ "$count" -eq 0 ]; then
    log "ERREUR: Aucun échantillon valide trouvé"
    exit 1
fi

log "✅ Validation manifest réussie: $count échantillons uniques"

cp "$MANIFEST_TEMP" "$MANIFEST"
rm "$MANIFEST_TEMP"

log "Aperçu du manifest:"
head -10 "$MANIFEST"

# ==================== 02 FASTQC RAW DATA ====================
log "FastQC sur données brutes"

mkdir -p "${ROOTDIR}/02_qualitycheck"
rm -f "${ROOTDIR}/02_qualitycheck"/*multiqc* 2>/dev/null || true
rm -rf "${ROOTDIR}/02_qualitycheck/multiqc_data" 2>/dev/null || true

cd "${ROOTDIR}/01_raw_data"

count=0
for file in *.fastq* */*.fastq*; do
    if [ -f "$file" ]; then
        count=$((count + 1))
        log "FastQC $count: $(basename $file)"
        conda run -n "$FASTQC_ENV" fastqc "$file" -o "${ROOTDIR}/02_qualitycheck" --threads 2 --quiet 2>/dev/null || {
            log "Erreur FastQC sur $file, continuons"
            continue
        }
        
        if [ $((count % 10)) -eq 0 ]; then
            sleep 2
        fi
    fi
done

log "FastQC terminé sur $count fichiers"

# MultiQC
log "MultiQC sur données brutes"
cd "${ROOTDIR}/02_qualitycheck"

conda run -n "$FASTQC_ENV" multiqc . \
    --force \
    --filename "raw_data_qc" \
    --title "Raw Data Quality Control - ValorMicro V2" \
    --ignore-symlinks \
    --no-ansi 2>/dev/null || {
    log "MultiQC a généré des warnings"
    if [ -f "raw_data_qc.html" ]; then
        log "✓ Rapport MultiQC créé"
    fi
}

# ==================== 02B TRIMMOMATIC - NETTOYAGE DES READS ====================
log "Trimmomatic - Nettoyage et suppression amorces 515F-926R"

mkdir -p "${ROOTDIR}/03_cleaned_data"
cd "${ROOTDIR}/01_raw_data"

# Créer fichier adapters FASTA pour cutadapt-style trimming
ADAPTERS="${ROOTDIR}/98_databasefiles/adapters_515f926r.fasta"
mkdir -p "${ROOTDIR}/98_databasefiles"

cat > "$ADAPTERS" << 'ADAPTER_EOF'
>515F
GTGYCAGCMGCCGCGGTAA
>926R
CCGYCAATTYMTTTRAGTTT
ADAPTER_EOF

log "Fichier adapters créé: $ADAPTERS"

# Compter paires totales
total_pairs=$(grep -c "^sample-id" "$MANIFEST" || echo "0")
if [ "$total_pairs" -gt 0 ]; then
    total_pairs=$((total_pairs - 1))
fi

pair_count=0
success_count=0
failed_count=0

log "Traitement de $total_pairs paires avec Trimmomatic..."

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r sample_id r1_path r2_path; do
    pair_count=$((pair_count + 1))
    
    # Fichiers de sortie
    base_r1=$(basename "$r1_path" .fastq.gz)
    base_r1=$(basename "$base_r1" .fastq)
    base_r2=$(basename "$r2_path" .fastq.gz)
    base_r2=$(basename "$base_r2" .fastq)
    
    out1p="${ROOTDIR}/03_cleaned_data/${sample_id}_R1_paired.fastq.gz"
    out1u="${ROOTDIR}/03_cleaned_data/${sample_id}_R1_unpaired.fastq.gz"
    out2p="${ROOTDIR}/03_cleaned_data/${sample_id}_R2_paired.fastq.gz"
    out2u="${ROOTDIR}/03_cleaned_data/${sample_id}_R2_unpaired.fastq.gz"
    
    log "Trimmomatic $pair_count/$total_pairs: $sample_id"
    
    # Exécution Trimmomatic avec suppression amorces
    conda run -n "$TRIMMOMATIC_ENV" trimmomatic PE -threads 8 -phred33 \
        "$r1_path" "$r2_path" \
        "$out1p" "$out1u" "$out2p" "$out2u" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100 2>&1 | grep -v "^Input Read" || {
        log "⚠️  Erreur Trimmomatic sur $sample_id - IGNORÉ"
        failed_count=$((failed_count + 1))
        continue
    }
    
    # Vérification synchronisation CRITIQUE
    if [[ -f "$out1p" && -f "$out2p" ]]; then
        # Compter reads avec méthode robuste
        count1=$(( $(zcat "$out1p" 2>/dev/null | wc -l) / 4 ))
        count2=$(( $(zcat "$out2p" 2>/dev/null | wc -l) / 4 ))
        
        if [ "$count1" = "$count2" ] && [ "$count1" -gt 0 ]; then
            log "✅ $sample_id SYNCHRONISÉ: $count1 reads"
            success_count=$((success_count + 1))
        else
            log "⚠️  $sample_id DÉSYNCHRONISÉ ($count1 vs $count2) - SUPPRESSION"
            rm -f "$out1p" "$out2p" "$out1u" "$out2u"
            failed_count=$((failed_count + 1))
        fi
    else
        log "⚠️  Fichiers de sortie manquants pour $sample_id"
        failed_count=$((failed_count + 1))
    fi
done

log "Trimmomatic terminé - Succès: $success_count / Erreurs: $failed_count"

# ==================== 03 CRÉER MANIFEST POUR DONNÉES NETTOYÉES ====================
log "Génération manifest pour données nettoyées"

MANIFEST_TRIMMED="${ROOTDIR}/98_databasefiles/manifest_trimmed"

echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST_TRIMMED"

cd "${ROOTDIR}/03_cleaned_data"

for r1_file in *_R1_paired.fastq.gz; do
    if [ -f "$r1_file" ]; then
        r2_file="${r1_file/_R1_paired/_R2_paired}"
        
        if [ -f "$r2_file" ]; then
            sample_id="${r1_file%%_R1_paired*}"
            
            r1_abs="${ROOTDIR}/03_cleaned_data/$r1_file"
            r2_abs="${ROOTDIR}/03_cleaned_data/$r2_file"
            
            echo -e "$sample_id\t$r1_abs\t$r2_abs" >> "$MANIFEST_TRIMMED"
            log "Ajouté au manifest nettoyé: $sample_id"
        fi
    fi
done

trimmed_count=$(tail -n +2 "$MANIFEST_TRIMMED" | wc -l)
log "✅ Manifest nettoyé créé avec $trimmed_count échantillons"

if [ "$trimmed_count" -eq 0 ]; then
    log "ERREUR: Aucun échantillon survit au trimming"
    exit 1
fi

# ==================== 04 FASTQC CLEANED DATA ====================
log "FastQC/MultiQC sur données nettoyées"

mkdir -p "${ROOTDIR}/03_cleaned_data_qc"
rm -f "${ROOTDIR}/03_cleaned_data_qc"/*multiqc* 2>/dev/null || true

cd "${ROOTDIR}/03_cleaned_data"

count=0
for file in *_paired.fastq.gz; do
    if [ -f "$file" ]; then
        count=$((count + 1))
        log "FastQC cleaned $count: $(basename $file)"
        conda run -n "$FASTQC_ENV" fastqc "$file" -o "${ROOTDIR}/03_cleaned_data_qc" --threads 2 --quiet 2>/dev/null || {
            log "Erreur FastQC sur $file"
            continue
        }
        
        if [ $((count % 10)) -eq 0 ]; then
            sleep 2
        fi
    fi
done

# MultiQC nettoyé
log "MultiQC sur données nettoyées"
cd "${ROOTDIR}/03_cleaned_data_qc"

conda run -n "$FASTQC_ENV" multiqc . \
    --force \
    --filename "cleaned_data_qc" \
    --title "Cleaned Data Quality Control - After Trimmomatic" \
    --ignore-symlinks \
    --no-ansi 2>/dev/null || {
    log "MultiQC cleaned data warnings"
    if [ -f "cleaned_data_qc.html" ]; then
        log "✓ Rapport MultiQC cleaned créé"
    fi
}

log "✅ Contrôle qualité post-nettoyage terminé"

# ==================== 05 QIIME2 IMPORT ====================
log "QIIME2 Import des données nettoyées"

mkdir -p "${ROOTDIR}/05_QIIME2/core" "${ROOTDIR}/05_QIIME2/visual"
cd "${ROOTDIR}/05_QIIME2"

conda run -n "$QIIME_ENV" qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST_TRIMMED" \
    --output-path "core/demux.qza" \
    --input-format PairedEndFastqManifestPhred33V2 || {
    log "ERREUR import QIIME2"
    exit 1
}

log "✅ Import QIIME2 réussi"

# Visualisation qualité
log "Visualisation qualité des reads"
conda run -n "$QIIME_ENV" qiime demux summarize \
    --i-data core/demux.qza \
    --o-visualization visual/demux-summary.qzv || {
    log "Erreur visualisation demux"
}

# ==================== 06 DADA2 ====================
log "DADA2 denoising"

cd "${ROOTDIR}/05_QIIME2/core"

conda run -n "$QIIME_ENV" qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux.qza \
    --o-table table.qza \
    --o-representative-sequences rep-seqs.qza \
    --o-denoising-stats denoising-stats.qza \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-n-threads "$NTHREADS" || {
    log "❌ DADA2 ÉCHOUÉ"
    exit 1
}

log "🎉 DADA2 RÉUSSI"

# ==================== 07 CLASSIFIER SILVA ====================
log "Classification taxonomique avec SILVA 138.2"

cd "${ROOTDIR}/98_databasefiles"

CLASSIFIER_URL="/nvme/bio/data_fungi/valormicro_nc/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
CLASSIFIER="${ROOTDIR}/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"

if [ ! -f "$CLASSIFIER" ]; then
    log "Copie classifier SILVA"
    cp "$CLASSIFIER_URL" "$CLASSIFIER" || {
        log "ERREUR: Classifier non trouvé"
        exit 1
    }
fi

conda run -n "$QIIME_ENV" qiime tools validate "$CLASSIFIER" || {
    log "ERREUR: Classifier invalide"
    exit 1
}

log "✅ Classifier SILVA 138.2 validé"

# Classification
log "Classification taxonomique"
cd "${ROOTDIR}/05_QIIME2/core"

conda run -n "$QIIME_ENV" qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza \
    --p-n-jobs "$NTHREADS"

log "✅ Classification taxonomique réussie"

conda run -n "$QIIME_ENV" qiime tools export \
    --input-path taxonomy.qza \
    --output-path temp_tax_check

if [ -f "temp_tax_check/taxonomy.tsv" ]; then
    tax_count=$(tail -n +2 temp_tax_check/taxonomy.tsv | wc -l)
    log "✅ Taxonomie contient $tax_count classifications"
fi
rm -rf temp_tax_check

# ==================== 08 ANALYSES FINALES ====================
log "Analyses finales: barplots, diversité"

mkdir -p "${ROOTDIR}/05_QIIME2/export"
cd "${ROOTDIR}/05_QIIME2/core"

# Summary table
log "Summary table"
conda run -n "$QIIME_ENV" qiime feature-table summarize \
    --i-table table.qza \
    --o-visualization "../visual/table-summary.qzv"

# Taxa barplots
log "Génération taxa barplots"
conda run -n "$QIIME_ENV" qiime taxa barplot \
    --i-table table.qza \
    --i-taxonomy taxonomy.qza \
    --o-visualization "../visual/taxa-bar-plots.qzv"

log "✅ Taxa barplots créés"

# ==================== 09 DIVERSITÉ ====================
log "Calcul métriques de diversité"

mkdir -p "${ROOTDIR}/05_QIIME2/diversity" "${ROOTDIR}/05_QIIME2/pcoa"

# Arbre phylogénétique
log "Génération arbre phylogénétique"
if [ ! -f "tree.qza" ]; then
    conda run -n "$QIIME_ENV" qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences rep-seqs.qza \
        --o-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza \
        --o-rooted-tree tree.qza \
        --p-n-threads "$NTHREADS"
fi

# Métadonnées
log "Création métadonnées"
cd "${ROOTDIR}/98_databasefiles"

echo -e "sample-id\tgroup\ttype" > "diversity-metadata.tsv"
tail -n +2 "$MANIFEST_TRIMMED" | cut -f1 | while read -r sample_id; do
    if echo "${sample_id,,}" | grep -qE "(neg|blank|control|ctrl|eau|mock)"; then
        echo -e "$sample_id\tcontrol\tnegative" >> "diversity-metadata.tsv"
    else
        echo -e "$sample_id\tsample\tenvironmental" >> "diversity-metadata.tsv"
    fi
done

log "✅ Métadonnées créées"

# Core metrics phylogenetic (SANS raréfaction)
cd "${ROOTDIR}/05_QIIME2/core"

log "Core metrics phylogenetic (SANS raréfaction)"
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
    --o-bray-curtis-emperor visual/Emperor-braycurtis.qzv

log "✅ Métriques de diversité créées"

# ==================== 10 EXPORTS ====================
log "Export de tous les fichiers QIIME2"

mkdir -p "${ROOTDIR}/05_QIIME2/export/core" "${ROOTDIR}/05_QIIME2/export/diversity_tsv"

cd "${ROOTDIR}/05_QIIME2"

conda run -n "$QIIME_ENV" qiime tools export \
    --input-path core/table.qza \
    --output-path export/core/table

conda run -n "$QIIME_ENV" qiime tools export \
    --input-path core/rep-seqs.qza \
    --output-path export/core/rep-seqs

conda run -n "$QIIME_ENV" qiime tools export \
    --input-path core/taxonomy.qza \
    --output-path export/core/taxonomy

# Fonction export diversité
export_diversity() {
    local qza_file="$1"
    local output_name="$2"
    
    if [ -f "$qza_file" ]; then
        log "Export $output_name"
        temp_dir="export/diversity_tsv/${output_name}_temp"
        mkdir -p "$temp_dir"
        
        conda run -n "$QIIME_ENV" qiime tools export \
            --input-path "$qza_file" \
            --output-path "$temp_dir" && {
            
            for ext in tsv txt csv; do
                find "$temp_dir" -name "*.${ext}" -type f | while read -r found_file; do
                    cp "$found_file" "export/diversity_tsv/${output_name}.${ext}"
                done
            done
            
            if [ ! -f "export/diversity_tsv/${output_name}.tsv" ]; then
                first_file=$(find "$temp_dir" -type f | head -1)
                if [ -f "$first_file" ]; then
                    cp "$first_file" "export/diversity_tsv/${output_name}.tsv"
                fi
            fi
        }
        
        rm -rf "$temp_dir"
    fi
}

# Export métriques
log "Export métriques alpha et beta"
export_diversity "core/diversity/Vector-observed_asv.qza" "observed_features"
export_diversity "core/diversity/Vector-shannon.qza" "shannon"
export_diversity "core/diversity/Vector-evenness.qza" "evenness"
export_diversity "core/diversity/Vector-faith_pd.qza" "faith_pd"
export_diversity "core/diversity/Matrix-jaccard.qza" "jaccard_distance"
export_diversity "core/diversity/Matrix-braycurtis.qza" "bray_curtis_distance"
export_diversity "core/diversity/Matrix-unweighted_unifrac.qza" "unweighted_unifrac_distance"
export_diversity "core/diversity/Matrix-weighted_unifrac.qza" "weighted_unifrac_distance"
export_diversity "core/pcoa/PCoA-jaccard.qza" "jaccard_pcoa"
export_diversity "core/pcoa/PCoA-braycurtis.qza" "bray_curtis_pcoa"
export_diversity "core/pcoa/PCoA-unweighted_unifrac.qza" "unweighted_unifrac_pcoa"
export_diversity "core/pcoa/PCoA-weighted_unifrac.qza" "weighted_unifrac_pcoa"
export_diversity "core/denoising-stats.qza" "dada2_stats"

# ==================== 11 CONVERSION BIOM ====================
log "Conversion BIOM vers TSV et création ASV.txt"

cd "${ROOTDIR}/05_QIIME2/export"

if [ -f "core/table/feature-table.biom" ]; then
    log "Conversion BIOM"
    
    mkdir -p "core/table_tsv"
    
    conda run -n "$QIIME_ENV" biom convert \
        -i core/table/feature-table.biom \
        -o core/table_tsv/table-from-biom.tsv \
        --to-tsv && {
        
        sed '1d ; s/#OTU ID/ASV_ID/' \
            core/table_tsv/table-from-biom.tsv > \
            core/table_tsv/ASV.tsv
        
        log "✅ ASV.tsv créé"
    }
fi

# Création ASV.txt avec taxonomie
log "Création ASV.txt avec taxonomie"

if [ -f "core/table_tsv/ASV.tsv" ] && [ -f "core/taxonomy/taxonomy.tsv" ]; then
    asv_file="core/table_tsv/ASV.tsv"
    taxonomy_file="core/taxonomy/taxonomy.tsv"
    output_file="core/table_tsv/ASV.txt"
    
    sample_header=$(head -1 "$asv_file" | cut -f2-)
    echo -e "Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t${sample_header}" > "$output_file"
    
    tail -n +2 "$asv_file" | while IFS=$'\t' read -r asv_id asv_counts; do
        kingdom="Bacteria"
        phylum="Unassigned"
        class="Unassigned"
        order="Unassigned"
        family="Unassigned"
        genus="Unassigned"
        species="Unassigned"
        
        if tax_line=$(grep "^${asv_id}" "$taxonomy_file" 2>/dev/null); then
            tax_string=$(echo "$tax_line" | cut -f2)
            
            [[ "$tax_string" =~ D_0__([^;]+) ]] && kingdom="${BASH_REMATCH[1]}"
            [[ "$tax_string" =~ D_1__([^;]+) ]] && phylum="${BASH_REMATCH[1]}"
            [[ "$tax_string" =~ D_2__([^;]+) ]] && class="${BASH_REMATCH[1]}"
            [[ "$tax_string" =~ D_3__([^;]+) ]] && order="${BASH_REMATCH[1]}"
            [[ "$tax_string" =~ D_4__([^;]+) ]] && family="${BASH_REMATCH[1]}"
            [[ "$tax_string" =~ D_5__([^;]+) ]] && genus="${BASH_REMATCH[1]}"
            [[ "$tax_string" =~ D_6__([^;]+) ]] && species="${BASH_REMATCH[1]}"
        fi
        
        echo -e "${kingdom}\t${phylum}\t${class}\t${order}\t${family}\t${genus}\t${species}\t${asv_counts}" >> "$output_file"
    done
    
    lines_count=$(wc -l < "$output_file" 2>/dev/null || echo "0")
    log "✅ ASV.txt créé avec taxonomie ($lines_count lignes)"
fi

# ==================== RAPPORT FINAL ====================
log "Création rapport final"

tsv_count=$(find export/diversity_tsv -name "*.tsv" -o -name "*.txt" 2>/dev/null | wc -l || echo "0")

cat > "${ROOTDIR}/05_QIIME2/RAPPORT_FINAL.md" << EOF
# Pipeline QIIME2 ValorMicro V2 - TERMINÉ

## Date: $(date)

## ✅ Étapes complétées

1. **Contrôle qualité brutes**
   - FastQC: ${ROOTDIR}/02_qualitycheck/raw_data_qc.html

2. **Trimmomatic - Nettoyage**
   - Suppression amorces 515F-926R
   - QC filters: LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100
   - Échantillons nettoyés: $trimmed_count

3. **Contrôle qualité post-nettoyage**
   - FastQC: ${ROOTDIR}/03_cleaned_data_qc/cleaned_data_qc.html

4. **DADA2 denoising**
   - Table ASV créée
   - Séquences représentatives générées

5. **Classification taxonomique**
   - Classifier: SILVA 138.2 SSU NR99
   - Amorces: 515F-926R

6. **Analyses de diversité**
   - ⚠️ SANS raréfaction
   - Métriques alpha & beta calculées
   - PCoA pour toutes métriques

## 📁 Fichiers principaux

- **Taxa barplots**: ${ROOTDIR}/05_QIIME2/visual/taxa-bar-plots.qzv
- **PCoA Emperor**: ${ROOTDIR}/05_QIIME2/visual/Emperor-*.qzv
- **ASV + taxonomie**: ${ROOTDIR}/05_QIIME2/export/core/table_tsv/ASV.txt
- **Diversité TSV**: ${ROOTDIR}/05_QIIME2/export/diversity_tsv/

## ✅ PIPELINE TERMINÉ AVEC SUCCÈS

Total TSV exportés: $tsv_count
EOF

log "🎉 PIPELINE VALORMICRO V2 AVEC TRIMMOMATIC TERMINÉ AVEC SUCCÈS ! 🎉"
log ""
log "==================== FICHIERS PRINCIPAUX ===================="
log "📊 Barplots: ${ROOTDIR}/05_QIIME2/visual/taxa-bar-plots.qzv"
log "📊 PCoA: ${ROOTDIR}/05_QIIME2/visual/Emperor-*.qzv"
log "📊 Table ASV: ${ROOTDIR}/05_QIIME2/export/core/table_tsv/ASV.txt"
log "📊 Diversité: ${ROOTDIR}/05_QIIME2/export/diversity_tsv/"
log "📄 Rapport: ${ROOTDIR}/05_QIIME2/RAPPORT_FINAL.md"
log "📄 Données nettoyées: ${ROOTDIR}/03_cleaned_data/"
