# SEM-based gene set analysis

Gene Set Analysis (GSA) via self-contained test for group effect on
signaling (directed) pathways based on SEM. The core of the methodology
is implemented in the RICF algorithm of
[`SEMrun()`](https://grassiMario.github.io/SEMgraph/reference/SEMrun.md),
recovering from RICF output node-specific group effect p-values, and
Brown’s combined permutation p-values of node activation and inhibition.

## Usage

``` r
SEMgsa(g = list(), data, group, method = "BH", alpha = 0.05, n_rep = 1000, ...)
```

## Arguments

- g:

  A list of pathways to be tested.

- data:

  A matrix or data.frame. Rows correspond to subjects, and columns to
  graph nodes (variables).

- group:

  A binary vector. This vector must be as long as the number of
  subjects. Each vector element must be 1 for cases and 0 for control
  subjects.

- method:

  Multiple testing correction method. One of the values available in
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). By default,
  method is set to "BH" (i.e., Benjamini-Hochberg correction).

- alpha:

  Gene set test significance level (default = 0.05).

- n_rep:

  Number of randomization replicates (default = 1000).

- ...:

  Currently ignored.

## Value

A list of 2 objects:

1.  "gsa", A data.frame reporting the following information for each
    pathway in the input list:

    - "No.nodes", pathway size (number of nodes);

    - "No.DEGs", number of differential espression genes (DEGs) within
      the pathway, after multiple test correction with one of the
      methods available in
      [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html);

    - "pert", pathway perturbation status (see details);

    - "pNA", Brown's combined P-value of pathway node activation;

    - "pNI", Brown's combined P-value of pathway node inhibition;

    - "PVAL", Bonferroni combined P-value of pNA, and pNI; i.e., 2\*
      min(pNA, PNI);

    - "ADJP", Adjusted Bonferroni P-value of pathway perturbation; i.e.,
      min(No.pathways \* PVAL; 1).

2.  "DEG", a list with DEGs names per pathways.

## Details

For gaining more biological insights into the functional roles of
pre-defined subsets of genes, node perturbation obtained from RICF
fitting has been combined with up- or down-regulation of genes from a
reference interactome to obtain overall pathway perturbation as follows:

- The node perturbation is defined as activated when the minimum among
  the p-values is positive; if negative, the status is inhibited.

- Up- or down- regulation of genes is computed from the weighted
  adjacency matrix of each pathway as column sum of weights(-1,0,1) over
  each source node. If the overall sum of node weights is below 1, the
  pathway is flagged as down-regulated, otherwise as up-regulated.

- The combination between these two quantities allows to define the
  direction (up or down) of gene perturbation. Up- or down regulated
  gene status, associated with node inhibition, indicates a decrease in
  activation (or increase in inhibition) in cases with respect to
  control group. Conversely, up- or down regulated gene status,
  associated with node activation, indicates an increase in activation
  (or decrease in inhibition) in cases with respect to control group.

## References

Grassi, M., Tarantino, B. (2022). SEMgsa: topology-based pathway
enrichment analysis with structural equation models. BMC Bioinformatics,
17 Aug, 23, 344. \<https://doi.org/10.1186/s12859-022-04884-8\>

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
# \dontrun{

# Nonparanormal(npn) transformation
als.npn <- transformData(alsData$exprs)$data
#> Conducting the nonparanormal transformation via shrunkun ECDF...done.

# Selection of FTD-ALS pathways from KEGG pathways

paths.name <- c("MAPK signaling pathway",
                "Protein processing in endoplasmic reticulum",
                "Endocytosis",
                "Wnt signaling pathway",
                "Neurotrophin signaling pathway",
                "Amyotrophic lateral sclerosis")

j <- which(names(kegg.pathways) %in% paths.name)

GSA <- SEMgsa(kegg.pathways[j], als.npn, alsData$group,
              method = "bonferroni", alpha = 0.05,
              n_rep = 1000)
#> k = 1 MAPK signaling pathway 
#> k = 2 Protein processing in endoplasmic reticulum 
#> k = 3 Endocytosis 
#> k = 4 Wnt signaling pathway 
#> k = 5 Neurotrophin signaling pathway 
#> k = 6 Amyotrophic lateral sclerosis 
GSA$gsa
#>                                             No.nodes No.DEGs     pert
#> Neurotrophin signaling pathway                   119      33 down act
#> Amyotrophic lateral sclerosis                    364      66   up inh
#> Endocytosis                                      252      59     <NA>
#> MAPK signaling pathway                           294      47   up act
#> Protein processing in endoplasmic reticulum      171      51   up act
#> Wnt signaling pathway                            166      34     <NA>
#>                                                      pNa          pNi
#> Neurotrophin signaling pathway              4.840572e-14 2.720046e-14
#> Amyotrophic lateral sclerosis               2.647549e-12 1.273426e-13
#> Endocytosis                                 3.013589e-11 2.664535e-12
#> MAPK signaling pathway                      6.145084e-12 4.988220e-09
#> Protein processing in endoplasmic reticulum 7.039591e-12 6.101954e-10
#> Wnt signaling pathway                       1.427281e-08 2.227794e-07
#>                                                     PVAL         ADJP
#> Neurotrophin signaling pathway              5.440093e-14 3.264056e-13
#> Amyotrophic lateral sclerosis               2.546852e-13 1.528111e-12
#> Endocytosis                                 5.329071e-12 3.197442e-11
#> MAPK signaling pathway                      1.229017e-11 7.374101e-11
#> Protein processing in endoplasmic reticulum 1.407918e-11 8.447509e-11
#> Wnt signaling pathway                       2.854561e-08 1.712737e-07
GSA$DEG
#> $`Neurotrophin signaling pathway`
#>  [1] "207"   "10000" "91860" "814"   "1399"  "2309"  "10818" "11213" "3667" 
#> [10] "3845"  "5604"  "4215"  "4217"  "5602"  "5600"  "5599"  "5601"  "4793" 
#> [19] "4893"  "4916"  "5296"  "5335"  "5336"  "5664"  "5781"  "5879"  "5906" 
#> [28] "387"   "8767"  "9252"  "6272"  "6654"  "7531" 
#> 
#> $`Amyotrophic lateral sclerosis`
#>  [1] "55860"  "10189"  "55626"  "9973"   "25978"  "1329"   "1337"   "1340"  
#>  [9] "1346"   "10540"  "51164"  "25981"  "146754" "27019"  "83544"  "10126" 
#> [17] "2081"   "2733"   "2882"   "10013"  "3178"   "3181"   "220988" "3309"  
#> [25] "3710"   "3799"   "3800"   "64837"  "89953"  "84557"  "81631"  "5608"  
#> [33] "5600"   "2475"   "55706"  "56901"  "4747"   "4741"   "9542"   "23511" 
#> [41] "23225"  "79902"  "4928"   "5532"   "5534"   "5688"   "5689"   "5718"  
#> [49] "5879"   "5903"   "6390"   "6391"   "6392"   "6396"   "81929"  "23064" 
#> [57] "10280"  "140775" "23435"  "84790"  "10382"  "84617"  "29979"  "29978" 
#> [65] "29796"  "55255" 
#> 
#> $Endocytosis
#>  [1] "116983" "10097"  "10096"  "57180"  "163"    "64411"  "55738"  "84364" 
#>  [9] "10565"  "10092"  "829"    "867"    "25978"  "128866" "51510"  "1212"  
#> [17] "27128"  "26052"  "2060"   "8729"   "9815"   "9146"   "3133"   "3134"  
#> [25] "3561"   "83737"  "3799"   "10015"  "8394"   "5338"   "23550"  "10890" 
#> [33] "9230"   "22841"  "57403"  "5867"   "5869"   "7879"   "9135"   "387"   
#> [41] "56904"  "4087"   "6643"   "58533"  "23111"  "51324"  "10254"  "7251"  
#> [49] "9559"   "51160"  "55737"  "51534"  "23325"  "9897"   "8976"   "147179"
#> [57] "644150" "11059"  "118813"
#> 
#> $`MAPK signaling pathway`
#>  [1] "207"       "208"       "774"       "775"       "776"       "9254"     
#>  [7] "836"       "1398"      "1436"      "1946"      "8822"      "2885"     
#> [13] "3480"      "8517"      "8681"      "4254"      "4296"      "7786"     
#> [19] "4215"      "9064"      "8491"      "1432"      "5599"      "23162"    
#> [25] "5601"      "2122"      "8569"      "4763"      "4772"      "4775"     
#> [31] "56034"     "5228"      "100137049" "5494"      "5532"      "5534"     
#> [37] "5567"      "84867"     "5879"      "5881"      "5906"      "5921"     
#> [43] "5922"      "9252"      "6654"      "23118"     "57551"    
#> 
#> $`Protein processing in endoplasmic reticulum`
#>  [1] "4287"   "9532"   "578"    "1410"   "8454"   "51009"  "3301"   "54788" 
#>  [9] "5611"   "285126" "55741"  "80267"  "2081"   "30001"  "26270"  "3320"  
#> [17] "3326"   "7184"   "3309"   "10960"  "11253"  "5602"   "5599"   "5601"  
#> [25] "8720"   "51360"  "7841"   "4780"   "23645"  "9978"   "6185"   "51128" 
#> [33] "6396"   "25956"  "29927"  "7095"   "11231"  "6400"   "64374"  "6500"  
#> [41] "6745"   "6747"   "7321"   "7322"   "7326"   "29979"  "29978"  "165324"
#> [49] "23190"  "7466"   "55432" 
#> 
#> $`Wnt signaling pathway`
#>  [1] "10297" "816"   "894"   "896"   "57680" "1452"  "80319" "23002" "1856" 
#> [10] "23291" "7976"  "8549"  "6885"  "4609"  "4772"  "4776"  "64840" "5532" 
#> [19] "5534"  "5567"  "5879"  "5881"  "9978"  "387"   "9475"  "84870" "8607" 
#> [28] "59343" "6477"  "6907"  "7089"  "7091"  "57216" "80326"
#> 

# }
```
