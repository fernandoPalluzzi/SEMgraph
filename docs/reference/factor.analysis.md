# Factor analysis for high dimensional data

Wrapper for Factor Analysis with potentially high dimensional variables
implement in the "cate" R package (Author: Jingshu Wang \[aut\],
Qingyuan Zhao \[aut, cre\] Maintainer: Qingyuan Zhao
\<qz280@cam.ac.uk\>) that is optimized for the high dimensional problem
where the number of samples n is less than the number of variables p.

## Usage

``` r
factor.analysis(Y, r = 1, method = "pc")
```

## Arguments

- Y:

  data matrix, a n\*p matrix

- r:

  number of factors (default, r =1)

- method:

  algorithm to be used, "pc" (default) or "ml"

## Value

a list of objects

- Gamma:

  estimated factor loadings

- Z:

  estimated latent factors

- Sigma:

  estimated noise variance matrix

## Details

The two methods extracted from "cate" are quasi-maximum likelihood (ml),
and principal component analysis (pc). The ml is iteratively solved the
EM algorithm using the PCA solution as the initial value. See Bai and Li
(2012) for more details.

## References

Jushan Bai and Kunpeng Li (2012). Statistical Analysis of Factor Models
of High Dimension. The Annals of Statistics, 40 (1), 436-465
\<https://doi.org/10.1214/11-AOS966\>

Jingshu Wang and Qingyuan Zhao (2020). cate: High Dimensional Factor
Analysis and Confounder Adjusted Testing and Estimation. R package
version 1.1.1. \<https://CRAN.R-project.org/package=cate\>

## Examples

``` r
# Nonparanormal(npn) transformation
als.npn <- transformData(alsData$exprs)$data
#> Conducting the nonparanormal transformation via shrunkun ECDF...done.

## pc
pc<- factor.analysis(Y = als.npn, r = 2, method = "pc")
head(pc$Gamma)
#>            [,1]        [,2]
#> [1,]  0.3023514  0.60012228
#> [2,]  0.7577627  0.52654004
#> [3,] -0.7554944 -0.26095445
#> [4,]  0.1151171 -0.70688315
#> [5,]  0.4806646 -0.06877752
#> [6,]  0.1340605 -0.58860977
head(pc$Z)
#>             [,1]         [,2]
#> [1,] -0.33026001 -0.869972921
#> [2,] -0.68976787 -1.361599830
#> [3,]  0.67371880  0.007571711
#> [4,] -1.39466075 -0.688934677
#> [5,]  0.07287238 -1.076990819
#> [6,] -0.99396433 -0.488201691
head(pc$Sigma)
#>       207       208     10000       284       285       317 
#> 0.5421869 0.1423013 0.3548809 0.4808143 0.7579812 0.6293163 

## ml
ml <- factor.analysis(Y = als.npn, r = 2, method = "ml")
head(ml$Gamma)
#>             [,1]        [,2]
#> 207    0.2641432  0.55929046
#> 208    0.7728752  0.53838279
#> 10000 -0.7715566 -0.27324054
#> 284    0.1192708 -0.69677472
#> 285    0.4539599 -0.04553308
#> 317    0.1480935 -0.55334494
head(ml$Z)
#>            [,1]       [,2]
#> ALS2 -0.3086763 -0.7838108
#> ALS3 -0.5541185 -1.2612329
#> ALS4  0.7607658  0.1601732
#> ALS5 -1.3418033 -0.5918722
#> ALS6  0.1523818 -0.9958831
#> ALS7 -1.0637363 -0.3341780
head(ml$Sigma)
#>       207       208     10000       284       285       317 
#> 0.6111499 0.1040704 0.3215149 0.4952729 0.7850598 0.6664684 
```
