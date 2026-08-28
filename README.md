# valormicro_V2
Characterization of marine microbial resources for analysis and enhancement of New Caledonia's natural heritage - Project from Drs **Anton Véronique** (CNRS/IRD Nouméa - New Caledonia)

This project is a continuation and extension of the Valormicro project, which led to the publication of this article: 

*Stenger, P.-L., Majorel, C., Valette, L., Ihage, W., Jardin-Camps, M., Jourand, P., & Anton-Leberre, V. (2026). Spatial structuring dominates over seasonality in tropical coastal microbiomes: Insights from New Caledonia's Indo-Pacific lagoon. Journal of Environmental Quality, 55, e70215. https://doi.org/10.1002/jeq2.70215*

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
