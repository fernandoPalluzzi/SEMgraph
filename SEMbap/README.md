# SEMbap: Bow-free covariance search and data de-correlation

*Mario Grassi and Barbara Tarantino*

Correspondence to: barbara.tarantino@unipv.it

**R session**

R version 4.1.0 (2021-05-18)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] matrixcalc_1.0-6     logger_0.2.2         org.Hs.eg.db_3.13.0  AnnotationDbi_1.54.1
 [5] IRanges_2.26.0       S4Vectors_0.30.2     Biobase_2.52.0       BiocGenerics_0.38.0 
 [9] wesanderson_0.3.6    SEMdata_0.1.2        SEMgraph_1.1.4       lavaan_0.6-12       
[13] igraph_1.3.5         ggplot2_1.0.1        data.table_1.14.4    dplyr_1.0.10        
[17] mvtnorm_1.1-3        lrpsadmm_0.2.0      

loaded via a namespace (and not attached):
 [1] httr_1.4.4             bit64_4.0.5            splines_4.1.0          assertthat_0.2.1      
 [5] blob_1.2.3             GenomeInfoDbData_1.2.6 yaml_2.3.6             robustbase_0.95-0     
 [9] pbivnorm_0.6.0         pillar_1.8.1           RSQLite_2.2.18         lattice_0.20-45       
[13] glue_1.6.2             RcppEigen_0.3.3.9.2    digest_0.6.30          XVector_0.32.0        
[17] RColorBrewer_1.1-3     colorspace_2.0-3       Matrix_1.5-1           plyr_1.8.7            
[21] pcaPP_2.0-2            pkgconfig_2.0.3        zlibbioc_1.38.0        scales_1.2.1          
[25] glasso_1.11            RSpectra_0.16-1        huge_1.3.5             tibble_3.1.8          
[29] KEGGREST_1.32.0        mgcv_1.8-41            generics_0.1.3         cachem_1.0.6          
[33] cli_3.4.1              mnormt_2.1.1           proto_1.0.0            magrittr_2.0.3        
[37] crayon_1.5.2           memoise_2.0.1          fansi_1.0.3            nlme_3.1-160          
[41] MASS_7.3-58.1          tools_4.1.0            lifecycle_1.0.3        stringr_1.4.1         
[45] munsell_0.5.0          Biostrings_2.60.2      compiler_4.1.0         GenomeInfoDb_1.28.4   
[49] rlang_1.0.6            RCurl_1.98-1.9         rstudioapi_0.14        cvTools_0.3.2         
[53] bitops_1.0-7           boot_1.3-28            gtable_0.3.1           DBI_1.1.3             
[57] reshape2_1.4.4         R6_2.5.1               fastmap_1.1.0          bit_4.0.4             
[61] utf8_1.2.2             stringi_1.7.8          Rcpp_1.0.9             vctrs_0.5.0           
[65] png_0.1-7              DEoptimR_1.0-11        tidyselect_1.2.0


## SEMbap_code_data folder structure

This folder contains the following data and files that can be used to reproduce the analysis of the manuscript. The folder structure can be summarised as follows:


 main_sim.R
 An R script to run Data simulations

 main_real.R
 An R script to run breast cancer data analysis

 help.R
 An R script with functions to import for running the analysis

 trrust_rawdata.human.txt
 Transcription factor data

 brca.RData
 Raw BRCA data

./csv/: simulations’ results

results_dense_100.csv
results_dense_400.csv
results_sparse_100.csv
results_sparse_400.csv


## Installation note:

    1. Set folder SEMbap_code_data as the current working directory.
    2. Specify a smaller number of iterations by argument “seed_vec” to reduce computing time for data simulations. 
