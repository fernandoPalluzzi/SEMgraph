# SEMgsa: topology-based pathway enrichment analysis with Structural Equation Models

*Mario Grassi and Barbara Tarantino*

Correspondence to: barbara.tarantino@unipv.it

This folder contains the following data and  R codes that can be used to reproduce the analysis of the manuscript:

Grassi M, Tarantino B. **SEMgsa: topology-based pathway enrichment analysis with structural equation models**. BMC Bioinformatics, 2022 Aug 17; 23(1):344. https://doi.org/10.1186/s12859-022-04884-8

The folder structure can be summarised as follows:

./data/:
    
    GSE172114_series_matrix.txt.gz
    GSE172114_rsem_gene_count_matrix_TMM_69samples.csv
    Raw COVID-19 source data

./R/:
    
    help.R
    An R script with functions to import for running the analysis
	
	main_sim.R
    An R script to run Data simulations
    
    main_true.R
    An R script to run True data analysis

## Notes
1. Set folder SEMgsa as the current working directory.
2. Specify a smaller number of iterations through the "N" argument to reduce computing time for data simulations.

## Data source
The original data used in this research can be downloaded at [**GSE172114**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172114).
