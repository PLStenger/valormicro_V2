# valormicro_V2
Characterization of marine microbial resources for analysis and enhancement of New Caledonia's natural heritage - Project from Drs **Véronique Anton** (CNRS/IRD Nouméa - New Caledonia)

This project is a continuation and extension of the Valormicro project, which led to the publication of this article: 

*Stenger, P.-L., Majorel, C., Valette, L., Ihage, W., Jardin-Camps, M., Jourand, P., & Anton-Leberre, V. (2026). Spatial structuring dominates over seasonality in tropical coastal microbiomes: Insights from New Caledonia's Indo-Pacific lagoon. Journal of Environmental Quality, 55, e70215. https://doi.org/10.1002/jeq2.70215*

## Interactive map of Valormicro V2 project

➡️ [Open the interactive map of Valormicro V2 sites](https://plstenger.github.io/99_Map_Valormicro_V2_sites.html)

### Installing pipeline :

First, open your terminal. Then, run these two command lines :

    cd -place_in_your_local_computer
    git clone https://github.com/PLStenger/valormicro_V2.git

### Update the pipeline in local by :

    git pull

### Run by :   

    # For first 'exploratory results'
    time nohup bash pipeline.sh &> pipeline.out
    time nohup bash pipeline_scratch.sh &> pipeline_scratch.out

    # Run this for the final pipeline
    time nohup bash 01_pipeline_QIIME2_PE_valormicro_V2.sh &>  01_pipeline_QIIME2_PE_valormicro_V2.out
    time nohup bash 02_rarefaction.sh &> 02_rarefaction.out
    time nohup bash 03_taxonomy.sh &> 03_taxonomy.out
    time nohup bash 03_taxonomy_eukaryota.sh &> 03_taxonomy_eukaryota.out
    time nohup bash 04_core_biom.sh &> 04_core_biom.out

    # JEDI pipeline
    time nohup bash 50_jedi_valormicro.sh &> 50_jedi_valormicro.out

    ## /!\ if Nextflow or Java problems with JEDI pipeline, please run before:

    conda create -y -n nextflow-java -c conda-forge openjdk=21
    conda activate nextflow-java
    java -version
    conda activate nextflow-java
    export JAVA_HOME="$CONDA_PREFIX"
    export NXF_JAVA_HOME="$CONDA_PREFIX"
    export JAVA_CMD="$CONDA_PREFIX/bin/java"
    export PATH="$CONDA_PREFIX/bin:$PATH"
    echo "CONDA_PREFIX = $CONDA_PREFIX"
    echo "JAVA_CMD = $JAVA_CMD"
    "$JAVA_CMD" -version
    nextflow -version

    # If problem with singularity, please run:

    cat > /nvme/bio/data_fungi/valormicro_V2/03_JEDI_pipeline/singularity_timeout.config <<'EOF'
    singularity {
        enabled = true
        pullTimeout = '2h'
        cacheDir = '/nvme/bio/data_fungi/valormicro_V2/03_JEDI_pipeline/container_cache'
    }
    
    process {
        withName: 'NFCORE_AMPLISEQ:AMPLISEQ:FORMATTAXONOMY' {
            container = 'docker://biocontainers/biocontainers:v1.2.0_cv1'
        }
    }
    EOF
    
    export SINGULARITY_CACHEDIR=/nvme/bio/data_fungi/valormicro_V2/03_JEDI_pipeline/container_cache
    mkdir -p "$SINGULARITY_CACHEDIR"
    
    singularity pull \
      "$SINGULARITY_CACHEDIR/biocontainers-v1.2.0_cv1.sif" \
      "docker://biocontainers/biocontainers:v1.2.0_cv1"
