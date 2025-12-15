
# SEMtree: tree-based structure learning methods with Structural Equation Models

*Mario Grassi and Barbara Tarantino*

Correspondence to: barbara.tarantino@unipv.it

This folder contains the following data and  R codes that can be used to reproduce the analysis of the manuscript:

Grassi M, Tarantino B. **SEMtree: tree-based structure learning methods with structural equation models**. Bioinformatics, 2023 June 09; 39(6):btad377. https://doi.org/10.1093/bioinformatics/btad377

&nbsp;

The folder structure can be summarised as follows:
  
./data/:
    
    C0206750_disease_gda_summary.tsv
	GSE172114_rsem_gene_count_matrix_TMM_69samples.csv
    GA_result.Rdata
    GOres.txt
    GSE172114_rsem_gene_count_matrix_TMM_69samples.csv
    GSE172114_series_matrix.txt.gz
    ijmsv19p0402s2.xlsx
	cdat.tab
	rdat.tab
	xdat.tab
	source and tab covid19 data 

./R/:
    
    tree_sim.R
    An R script to run Data simulations
    
    tree_bench.R
    An R script to run True data analysis
    
    tree_comparison.R
    An R script to run Method comparison
    
    Help.R
    An R script with functions to import for running the analysis

## Notes
1. Set SEMtree as the current working directory.
2. To reduce computing time for data simulations, reduce the number of iterations through the argument “seed_vec”.

## Data source
The original data used in this research can be downloaded at [**GSE172114**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172114).
