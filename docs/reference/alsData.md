# Amyotrophic Lateral Sclerosis (ALS) dataset

Expression profiling through high-throughput sequencing (RNA-seq) of 139
ALS patients and 21 healthy controls (HCs), from Tam et al. (2019).

## Usage

``` r
alsData
```

## Format

alsData is a list of 4 objects:

1.  "graph", ALS graph as the largest connected component of the
    "Amyotrophic lateral sclerosis (ALS)" pathway from KEGG database;

2.  "exprs", a matrix of 160 rows (subjects) and 318 columns (genes)
    extracted from the original 17695. This subset includes genes from
    KEGG pathways, needed to run SEMgraph examples. Raw data from the
    GEO dataset GSE124439 (Tam et al., 2019) were pre-processed applying
    batch effect correction, using the sva R package (Leek et al.,
    2012), to remove data production center and brain area biases. Using
    multidimensional scaling-based clustering, ALS-specific and an
    HC-specific clusters were generated. Misclassified samples were
    blacklisted and removed from the current dataset;

3.  "group", a binary group vector of 139 ALS subjects (1) and 21
    healthy controls (0);

4.  "details", a data.frame reporting information about included and
    blacklisted samples.

## Source

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124439>

## References

Tam OH, Rozhkov NV, Shaw R, Kim D et al. (2019). Postmortem Cortex
Samples Identify Distinct Molecular Subtypes of ALS: Retrotransposon
Activation, Oxidative Stress, and Activated Glia. Cell Reports,
29(5):1164-1177.e5. \<https://doi.org/10.1016/j.celrep.2019.09.066\>

Jeffrey T. Leek, W. Evan Johnson, Hilary S. Parker, Andrew E. Jaffe, and
John D. Storey (2012). The sva package for removing batch effects and
other unwanted variation in high-throughput experiments. Bioinformatics.
Mar 15; 28(6): 882-883.
\<https://doi.org/10.1093/bioinformatics/bts034\>

## Examples

``` r
alsData$graph
#> IGRAPH 5d44052 DNW- 32 47 -- 
#> + attr: name (v/c), weight (e/n)
#> + edges from 5d44052 (vertex names):
#>  [1] 6647 ->10452 6647 ->84134 6647 ->596   6647 ->4747  6647 ->79139
#>  [6] 6647 ->5530  6647 ->5532  6647 ->5533  6647 ->5534  6647 ->5535 
#> [11] 54205->842   7124 ->7132  7124 ->7133  581  ->54205 572  ->54205
#> [16] 596  ->54205 598  ->54205 317  ->842   842  ->836   7132 ->1616 
#> [21] 7133 ->1616  1616 ->4217  4217 ->5606  4217 ->5608  5606 ->1432 
#> [26] 5606 ->5600  5606 ->5603  5606 ->6300  5608 ->1432  5608 ->5600 
#> [31] 5608 ->5603  5608 ->6300  1432 ->4747  1432 ->4741  1432 ->4744 
#> [36] 5600 ->4747  5600 ->4741  5600 ->4744  5603 ->4747  5603 ->4741 
#> + ... omitted several edges
dim(alsData$exprs)
#> [1] 160 318
table(alsData$group)
#> 
#>   0   1 
#>  21 139 
```
