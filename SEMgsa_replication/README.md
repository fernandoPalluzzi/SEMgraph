# SEMgsa: topology-based pathway enrichment analysis with Structural Equation Models

*by Mario Grassi and Barbara Tarantino*

This folder contains the following data and files that can be used to
reproduce the analysis of the manuscript. The folder structure can be
summarised as follows:

./Data/:
    
    GSE172114_series_matrix.txt.gz
    GSE172114_rsem_gene_count_matrix_TMM_69samples.csv
    Raw COVID-19 source data

./R/:
    
    Main_sim.R
    An R script to run Data simulations
    
    Main_true.R
    An R script to run True data analysis

    Help.R
    An R script with functions to import for running the analysis

## Notes
1. Set folder SEMgsa_code_data as the current working directory.
2. Specify a smaller number of iterations through the "N" argument to reduce computing time for data simulations.

## Data source
The original data used in this research can be downloaded at [**GSE172114**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172114).

## References

Mario Grassi, Barbara Tarantino. **SEMgsa: topology-based pathway enrichment analysis with structural equation models.** BMC Bioinformatics. 2022 Aug 17;23(1):344. doi: 10.1186/s12859-022-04884-8.
