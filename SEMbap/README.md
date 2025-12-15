# SEMbap: Bow-free covariance search and data de-correlation

*Mario Grassi and Barbara Tarantino*

Correspondence to: barbara.tarantino@unipv.it

This folder contains the following data and  R codes that can be used to reproduce the analysis of the manuscript:

Grassi M, Tarantino B. **SEMbap: Bow-free covariance search and data de-correlation**. PLoS Comput Biol, 2024 Sep 11; 20(9):e1012448. https://doi.org/10.1371/journal.pcbi.1012448

&nbsp;

The folder structure can be summarised as follows:
  
./data/:
 
    brca.RData
    Raw BRCA data
	
	trrust_rawdata.human.txt
    Transcription factor data

    ./csv/: 

    results_dense_100.csv
    results_dense_400.csv
    results_sparse_100.csv
    results_sparse_400.csv
    simulations’ results

./R/:

    help.R
    An R script with functions to import for running the analysis
	
	main_sim.R
    An R script to run data simulations

    main_real.R
    An R script to run breast cancer data analysis

## Installation note:

1. Set folder SEMbap as the current working directory.
2. Specify a smaller number of iterations by argument “seed_vec” to reduce computing time for data simulations. 
