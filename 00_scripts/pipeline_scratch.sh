#!/usr/bin/env bash

set -euo pipefail

# ==================== CONFIGURATION ====================
export ROOTDIR="/scratch_vol0/fungi/valormicro_V2"
export NTHREADS=16
export TMPDIR="${ROOTDIR}/tmp"

mkdir -p "$TMPDIR"

log() { echo -e "\n[$(date +'%F %T')] $*\n"; }

log "=== PIPELINE VALORMICRO V2 DÉMARRÉ ==="

# ---- CORRECTION TOKEN CONDA ----
log "Suppression token conda corrompu"
sudo rm -rf /home/fungi/.conda/ 2>/dev/null || rm -rf /home/fungi/.conda/ 2>/dev/null || true

set +u
source $(conda info --base)/etc/profile.d/conda.sh
export JAVA_HOME="${JAVA_HOME:-/usr/lib/jvm/default-java}"
set -u

# ==================== 01 GÉNÉRATION MANIFEST ====================
log "Génération manifest pour tous les échantillons"

cd "${ROOTDIR}"
mkdir -p 98_databasefiles

MANIFEST="${ROOTDIR}/98_databasefiles/manifest"

# Créer header
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST"

# Scanner fichiers R1/R2
cd "${ROOTDIR}/01_raw_data"

log "Scan des fichiers fastq dans: $(pwd)"

declare -A seen_ids
count=0

for r1_file in *R1*.fastq* *_R1.fastq*; do
    if [ -f "$r1_file" ]; then
        # Trouver R2 correspondant
        r2_file="${r1_file/R1/R2}"
        r2_file="${r2_file/_R1./_R2.}"
        
        if [ -f "$r2_file" ]; then
            # Vérifier taille des fichiers (>1KB)
            r1_size=$(stat -c%s "$r1_file" 2>/dev/null || echo "0")
            r2_size=$(stat -c%s "$r2_file" 2>/dev/null || echo "0")
            
            if [ "$r1_size" -gt 1000 ] && [ "$r2_size" -gt 1000 ]; then
                # Extraire sample-id
                base_name=$(basename "$r1_file")
                sample_id="${base_name%%_S[0-9]*}"  # Enlever _S## et tout ce qui suit
                sample_id="${sample_id%%_R1*}"      # Enlever _R1 et ce qui suit
                sample_id="${sample_id%.fastq*}"    # Enlever extensions
                
                # Remplacer caractères problématiques
                sample_id="${sample_id//[^a-zA-Z0-9._-]/_}"
                
                # Gérer doublons
                original_id="$sample_id"
                counter=1
                while [[ -n "${seen_ids[$sample_id]:-}" ]]; do
                    sample_id="${original_id}_${counter}"
                    counter=$((counter + 1))
                done
                seen_ids["$sample_id"]=1
                
                # Chemins absolus
                r1_abs="${ROOTDIR}/01_raw_data/$r1_file"
                r2_abs="${ROOTDIR}/01_raw_data/$r2_file"
                
                echo -e "$sample_id\t$r1_abs\t$r2_abs" >> "$MANIFEST"
                count=$((count + 1))
                log "Ajouté: $sample_id"
            fi
        fi
    fi
done

log "Manifest créé avec $count échantillons"

if [ "$count" -eq 0 ]; then
    log "ERREUR: Aucun échantillon trouvé"
    exit 1
fi

# Vérifier doublons
duplicates=$(cut -f1 "$MANIFEST" | sort | uniq -d)
if [ -n "$duplicates" ]; then
    log "ERREUR: Doublons détectés: $duplicates"
    exit 1
fi

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
        conda run -n fastqc fastqc "$file" -o "${ROOTDIR}/02_qualitycheck" --threads 2 --quiet || {
            log "Erreur FastQC sur $file, continuons"
            continue
        }
        
        # Pause tous les 10 fichiers
        if [ $((count % 10)) -eq 0 ]; then
            sleep 2
        fi
    fi
done

log "FastQC terminé sur $count fichiers"

# MultiQC
log "MultiQC sur données brutes"
cd "${ROOTDIR}/02_qualitycheck"

conda run -n multiqc multiqc . \
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

# ==================== 03 QIIME2 IMPORT ====================
log "QIIME2 Import des données brutes"

mkdir -p "${ROOTDIR}/05_QIIME2/core" "${ROOTDIR}/05_QIIME2/visual"
cd "${ROOTDIR}/05_QIIME2"

conda run -n qiime2-amplicon-2025.7 qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST" \
    --output-path "core/demux.qza" \
    --input-format PairedEndFastqManifestPhred33V2 || {
    log "ERREUR import QIIME2"
    exit 1
}

log "✅ Import QIIME2 réussi"

# Visualisation qualité
log "Visualisation qualité des reads"
conda run -n qiime2-amplicon-2025.7 qiime demux summarize \
    --i-data core/demux.qza \
    --o-visualization visual/demux-summary.qzv || {
    log "Erreur visualisation demux"
}

# ==================== 04 CUTADAPT - SUPPRESSION AMORCES ====================
log "Cutadapt - Suppression amorces 515F-926R"

# Amorces du script original (région V4-V5)
PRIMER_F="GTGYCAGCMGCCGCGGTAA"
PRIMER_R="CCGYCAATTYMTTTRAGTTT"

conda run -n qiime2-amplicon-2025.7 qiime cutadapt trim-paired \
    --i-demultiplexed-sequences core/demux.qza \
    --p-front-f "$PRIMER_F" \
    --p-front-r "$PRIMER_R" \
    --p-error-rate 0.1 \
    --p-overlap 3 \
    --p-match-read-wildcards \
    --p-match-adapter-wildcards \
    --p-discard-untrimmed \
    --o-trimmed-sequences core/demux-trimmed.qza \
    --verbose || {
    log "ERREUR Cutadapt"
    exit 1
}

log "✅ Cutadapt terminé"

# Visualisation après cutadapt
log "Visualisation après suppression amorces"
conda run -n qiime2-amplicon-2025.7 qiime demux summarize \
    --i-data core/demux-trimmed.qza \
    --o-visualization visual/demux-trimmed-summary.qzv || {
    log "Erreur visualisation trimmed"
}

# ==================== 05 FASTQC CLEANED DATA ====================
log "FastQC/MultiQC sur données après suppression amorces"

mkdir -p "${ROOTDIR}/03_cleaned_data_qc"
rm -f "${ROOTDIR}/03_cleaned_data_qc"/*multiqc* 2>/dev/null || true

# Export temporaire pour FastQC
log "Export temporaire pour FastQC"
mkdir -p "${ROOTDIR}/temp_export_cleaned"

conda run -n qiime2-amplicon-2025.7 qiime tools export \
    --input-path core/demux-trimmed.qza \
    --output-path "${ROOTDIR}/temp_export_cleaned" || {
    log "Erreur export pour FastQC"
}

# FastQC sur fichiers nettoyés
cd "${ROOTDIR}/temp_export_cleaned"
count=0
for file in *.fastq* */*.fastq*; do
    if [ -f "$file" ]; then
        count=$((count + 1))
        log "FastQC cleaned $count: $(basename $file)"
        conda run -n fastqc fastqc "$file" -o "${ROOTDIR}/03_cleaned_data_qc" --threads 2 --quiet || {
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

conda run -n multiqc multiqc . \
    --force \
    --filename "cleaned_data_qc" \
    --title "Cleaned Data Quality Control - After Primer Removal" \
    --ignore-symlinks \
    --no-ansi 2>/dev/null || {
    log "MultiQC cleaned data warnings"
    if [ -f "cleaned_data_qc.html" ]; then
        log "✓ Rapport MultiQC cleaned créé"
    fi
}

# Nettoyer export temporaire
rm -rf "${ROOTDIR}/temp_export_cleaned"

log "✅ Contrôle qualité post-nettoyage terminé"

# ==================== 06 DADA2 ====================
log "DADA2 denoising"

cd "${ROOTDIR}/05_QIIME2/core"

conda run -n qiime2-amplicon-2025.7 qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux-trimmed.qza \
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

# URL du classifier du script original
CLASSIFIER_URL="/scratch_vol0/fungi/valormicro_nc/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
CLASSIFIER="${ROOTDIR}/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"

if [ ! -f "$CLASSIFIER" ]; then
    log "Téléchargement classifier SILVA depuis l'adresse du script"
    wget -O "$CLASSIFIER" "$CLASSIFIER_URL" || {
        log "ERREUR: Téléchargement classifier échoué"
        log "Le classifier doit être disponible à: $CLASSIFIER_URL"
        exit 1
    }
fi

# Valider classifier
conda run -n qiime2-amplicon-2025.7 qiime tools validate "$CLASSIFIER" || {
    log "ERREUR: Classifier invalide"
    exit 1
}

log "✅ Classifier SILVA 138.2 validé"

# Classification
log "Classification taxonomique"
cd "${ROOTDIR}/05_QIIME2/core"

conda run -n qiime2-amplicon-2025.7 qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza \
    --p-n-jobs "$NTHREADS"

log "✅ Classification taxonomique SILVA 138.2 réussie"

# Vérification taxonomie
conda run -n qiime2-amplicon-2025.7 qiime tools export \
    --input-path taxonomy.qza \
    --output-path temp_tax_check

if [ -f "temp_tax_check/taxonomy.tsv" ]; then
    tax_count=$(tail -n +2 temp_tax_check/taxonomy.tsv | wc -l)
    log "✅ Taxonomie contient $tax_count classifications"
    log "Aperçu taxonomie:"
    head -5 temp_tax_check/taxonomy.tsv | column -t -s$'\t' || head -5 temp_tax_check/taxonomy.tsv
fi
rm -rf temp_tax_check

# ==================== 08 ANALYSES FINALES ====================
log "Analyses finales: barplots, PCA, diversité"

mkdir -p "${ROOTDIR}/05_QIIME2/export"
cd "${ROOTDIR}/05_QIIME2/core"

# Summary table
log "Summary table"
conda run -n qiime2-amplicon-2025.7 qiime feature-table summarize \
    --i-table table.qza \
    --o-visualization "../visual/table-summary.qzv"

# Export summary
conda run -n qiime2-amplicon-2025.7 qiime tools export \
    --input-path "../visual/table-summary.qzv" \
    --output-path "../visual/table-summary"

# Taxa barplots
log "Génération taxa barplots"
conda run -n qiime2-amplicon-2025.7 qiime taxa barplot \
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
    conda run -n qiime2-amplicon-2025.7 qiime phylogeny align-to-tree-mafft-fasttree \
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
tail -n +2 "$MANIFEST" | cut -f1 | while read -r sample_id; do
    if echo "${sample_id,,}" | grep -qE "(neg|blank|control|ctrl|eau|mock)"; then
        echo -e "$sample_id\tcontrol\tnegative" >> "diversity-metadata.tsv"
    else
        echo -e "$sample_id\tsample\tenvironmental" >> "diversity-metadata.tsv"
    fi
done

log "✅ Métadonnées créées"

# Core metrics phylogenetic (SANS raréfaction)
cd "${ROOTDIR}/05_QIIME2/core"

log "Core metrics phylogenetic"
conda run -n qiime2-amplicon-2025.7 qiime diversity core-metrics-phylogenetic \
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

mkdir -p "${ROOTDIR}/05_QIIME2/export/core" \
         "${ROOTDIR}/05_QIIME2/export/diversity_tsv"

cd "${ROOTDIR}/05_QIIME2"

# Export table principale
conda run -n qiime2-amplicon-2025.7 qiime tools export \
    --input-path core/table.qza \
    --output-path export/core/table

# Export séquences
conda run -n qiime2-amplicon-2025.7 qiime tools export \
    --input-path core/rep-seqs.qza \
    --output-path export/core/rep-seqs

# Export taxonomie
conda run -n qiime2-amplicon-2025.7 qiime tools export \
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
        
        conda run -n qiime2-amplicon-2025.7 qiime tools export \
            --input-path "$qza_file" \
            --output-path "$temp_dir" && {
            
            # Trouver fichiers TSV/TXT
            for ext in tsv txt csv; do
                find "$temp_dir" -name "*.${ext}" -type f | while read -r found_file; do
                    cp "$found_file" "export/diversity_tsv/${output_name}.${ext}"
                    log "✅ ${output_name}.${ext} créé"
                done
            done
            
            # Fallback premier fichier
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

# Export métriques alpha
log "Export métriques alpha"
export_diversity "core/diversity/Vector-observed_asv.qza" "observed_features"
export_diversity "core/diversity/Vector-shannon.qza" "shannon"
export_diversity "core/diversity/Vector-evenness.qza" "evenness"
export_diversity "core/diversity/Vector-faith_pd.qza" "faith_pd"

# Export matrices distance
log "Export matrices distance"
export_diversity "core/diversity/Matrix-jaccard.qza" "jaccard_distance"
export_diversity "core/diversity/Matrix-braycurtis.qza" "bray_curtis_distance"
export_diversity "core/diversity/Matrix-unweighted_unifrac.qza" "unweighted_unifrac_distance"
export_diversity "core/diversity/Matrix-weighted_unifrac.qza" "weighted_unifrac_distance"

# Export PCoA
log "Export PCoA"
export_diversity "core/pcoa/PCoA-jaccard.qza" "jaccard_pcoa"
export_diversity "core/pcoa/PCoA-braycurtis.qza" "bray_curtis_pcoa"
export_diversity "core/pcoa/PCoA-unweighted_unifrac.qza" "unweighted_unifrac_pcoa"
export_diversity "core/pcoa/PCoA-weighted_unifrac.qza" "weighted_unifrac_pcoa"

# Export stats DADA2
export_diversity "core/denoising-stats.qza" "dada2_stats"

# ==================== 11 CONVERSION BIOM ====================
log "Conversion BIOM vers TSV et création ASV.txt"

cd "${ROOTDIR}/05_QIIME2/export"

# Conversion table
if [ -f "core/table/feature-table.biom" ]; then
    log "Conversion BIOM"
    
    mkdir -p "core/table_tsv"
    
    conda run -n qiime2-amplicon-2025.7 biom convert \
        -i core/table/feature-table.biom \
        -o core/table_tsv/table-from-biom.tsv \
        --to-tsv && {
        
        # Modifier header
        sed '1d ; s/#OTU ID/ASV_ID/' \
            core/table_tsv/table-from-biom.tsv > \
            core/table_tsv/ASV.tsv
        
        log "✅ ASV.tsv créé ($(wc -l < core/table_tsv/ASV.tsv) lignes)"
    }
fi

# Création ASV.txt avec taxonomie
log "Création ASV.txt avec taxonomie"

if [ -f "core/table_tsv/ASV.tsv" ] && [ -f "core/taxonomy/taxonomy.tsv" ]; then
    asv_file="core/table_tsv/ASV.tsv"
    taxonomy_file="core/taxonomy/taxonomy.tsv"
    output_file="core/table_tsv/ASV.txt"
    
    # Header
    sample_header=$(head -1 "$asv_file" | cut -f2-)
    echo -e "Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t${sample_header}" > "$output_file"
    
    # Traitement ASVs
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
    
    log "Aperçu ASV.txt:"
    head -3 "$output_file" | column -t -s$'\t' 2>/dev/null || head -3 "$output_file"
fi

# ==================== RAPPORT FINAL ====================
log "Création rapport final"

tsv_count=$(find export/diversity_tsv -name "*.tsv" -o -name "*.txt" 2>/dev/null | wc -l || echo "0")

cat > "${ROOTDIR}/05_QIIME2/RAPPORT_FINAL.md" << EOF
# Pipeline QIIME2 ValorMicro V2 - TERMINÉ

## Date: $(date)

## ✅ Étapes complétées

1. **Contrôle qualité**
   - FastQC sur données brutes: ${ROOTDIR}/02_qualitycheck/raw_data_qc.html
   - FastQC sur données nettoyées: ${ROOTDIR}/03_cleaned_data_qc/cleaned_data_qc.html

2. **Import et nettoyage QIIME2**
   - Import: $count échantillons
   - Suppression amorces: 515F-926R (région V4-V5)

3. **DADA2 denoising**
   - Table ASV créée
   - Séquences représentatives générées

4. **Classification taxonomique**
   - Classifier: SILVA 138.2 SSU NR99
   - Amorces: 515F-926R
   - Classifications: $tax_count ASVs

5. **Analyses de diversité**
   - ⚠️ SANS raréfaction (toutes données conservées)
   - Métriques alpha: richesse, Shannon, équitabilité, Faith PD
   - Métriques beta: Jaccard, Bray-Curtis, UniFrac
   - PCoA calculées pour toutes métriques beta

## 📁 Fichiers de sortie principaux

### Visualisations
- **Taxa barplots**: ${ROOTDIR}/05_QIIME2/visual/taxa-bar-plots.qzv
- **PCoA Emperor plots**: ${ROOTDIR}/05_QIIME2/visual/Emperor-*.qzv
- **Summary table**: ${ROOTDIR}/05_QIIME2/visual/table-summary.qzv

### Tables
- **ASV avec taxonomie**: ${ROOTDIR}/05_QIIME2/export/core/table_tsv/ASV.txt
- **Taxonomie**: ${ROOTDIR}/05_QIIME2/export/core/taxonomy/taxonomy.tsv
- **Table BIOM**: ${ROOTDIR}/05_QIIME2/export/core/table/feature-table.biom

### Diversité (TSV)
- **Métriques alpha**: ${ROOTDIR}/05_QIIME2/export/diversity_tsv/
  - observed_features.tsv
  - shannon.tsv
  - evenness.tsv
  - faith_pd.tsv
- **Matrices distance**: 
  - jaccard_distance.tsv
  - bray_curtis_distance.tsv
  - unweighted_unifrac_distance.tsv
  - weighted_unifrac_distance.tsv
- **PCoA**:
  - jaccard_pcoa.tsv
  - bray_curtis_pcoa.tsv
  - unweighted_unifrac_pcoa.tsv
  - weighted_unifrac_pcoa.tsv

### Stats
- **DADA2**: ${ROOTDIR}/05_QIIME2/export/diversity_tsv/dada2_stats.tsv

## 📊 Utilisation des fichiers

### Pour barplots de composition
Utiliser: \`visual/taxa-bar-plots.qzv\` (ouvrir sur view.qiime2.org)
Ou: \`export/core/table_tsv/ASV.txt\` pour graphiques personnalisés

### Pour PCA (PCoA)
Utiliser les fichiers: \`export/diversity_tsv/*_pcoa.tsv\`
Coordonnées prêtes pour visualisation

### Pour valeurs de diversité
Utiliser: \`export/diversity_tsv/observed_features.tsv\`, \`shannon.tsv\`, etc.

## ⚠️ NOTES IMPORTANTES

**Absence de raréfaction:**
- Les analyses de diversité utilisent TOUTES les séquences
- Les métriques de richesse/diversité reflètent la profondeur réelle de séquençage
- Les résultats peuvent être influencés par les différences de couverture entre échantillons
- À adapter selon vos besoins biologiques (filtrages additionnels possibles en post-traitement)

## ✅ PIPELINE TERMINÉ AVEC SUCCÈS

Total fichiers TSV exportés: $tsv_count
EOF

log "🎉 PIPELINE VALORMICRO V2 TERMINÉ AVEC SUCCÈS ! 🎉"
log ""
log "==================== FICHIERS PRINCIPAUX ===================="
log "📊 Barplots: ${ROOTDIR}/05_QIIME2/visual/taxa-bar-plots.qzv"
log "📊 PCoA: ${ROOTDIR}/05_QIIME2/visual/Emperor-*.qzv"
log "📊 Table ASV: ${ROOTDIR}/05_QIIME2/export/core/table_tsv/ASV.txt"
log "📊 Diversité TSV: ${ROOTDIR}/05_QIIME2/export/diversity_tsv/"
log "📄 Rapport: ${ROOTDIR}/05_QIIME2/RAPPORT_FINAL.md"
log ""
log "Visualiser les .qzv sur: https://view.qiime2.org"
