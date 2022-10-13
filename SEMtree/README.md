
# SEMtree: tree-based structure learning methods with Structural Equation Models

*by Mario Grassi and Barbara Tarantino*

*R code corresponding author: **Barbara Tarantino** (barbara.tarantino@unipv.it)*

&nbsp;

This section contains R code and data for the SEMtree project reproducibility.

The folder structure can be summarised as follows:

./code/:
    
    tree_sim.R
    An R script to run Data simulations
    
    tree_bench.R
    An R script to run True data analysis
    
    tree_comparison.R
    An R script to run Method comparison
    
    Help.R
    An R script with functions to import for running the analysis
    
./data/:
    
    C0206750_disease_gda_summaryGSE172114_rsem_gene_count_matrix_TMM_69samples.csv
    GA_result.Rdata
    Gores.txt
    GSE172114_rsem_gene_count_matrix_TMM_69samples.csv
    GSE172114_series_matrix.txt.gz
    ijmsv19p0402s2.xlsx

./graph/:
    
    BioNet
    COSINE
    KEGGCovid19.graph
    pathfinder
    ST
    STcfa
    STcov
    Stone
    STr2z
    STr2zN
    STsem
    tree_base
    WalktrapGM
    WalktrapGMl

## Notes
1. Set folder SEMtree_code_data as the current working directory.
2. Specify a smaller number of iterations by argument “seed_vec” to reduce computing time for data simulations.

## Data source
The original data used in this research can be downloaded at [**GSE172114**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172114).

## References

...
