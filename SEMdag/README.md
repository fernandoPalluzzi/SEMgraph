# SEMdag: fast learning of Directed Acyclic Graphs via node or layer ordering

*Mario Grassi and Barbara Tarantino*

Correspondence to: barbara.tarantino@unipv.it

This folder contains the following data and  R codes that can be used to reproduce the analysis of the manuscript:

Grassi M, Tarantino B. **SEMdag: Fast learning of Directed Acyclic Graphs via node or layer ordering**. PLoS ONE. 2025 Jan 08; 20(1): e0317283. https://doi.org/10.1371/journal.pone.0317283

&nbsp;

The folder structure can be summarised as follows:
  
./data/:
      als/
	  brca/
	  covid19/
	  stemi/
 	  train and test dataset

./R/:
    direct_lingam_funcs.cpp
	C++ functions for DirectLINGAM method
	
    DLinGAM_p.R
	An R script with DirectLINGAM function
	
    help.R
	An R script with functions to import for running the analysis
    
	main.R
	An R script to run data analysis

## Installation note:

1. Set folder SEMdag as the current working directory.
2. Specify a smaller number of iterations by argument “seed_vec” to reduce computing time for data simulations. 
