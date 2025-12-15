# RICF model summary

Generate a summary for a RICF solver similar to lavaan-formatted summary

## Usage

``` r
# S3 method for class 'RICF'
summary(object, ...)
```

## Arguments

- object:

  A RICF fitted model object.

- ...:

  Currently ignored.

## Value

Shown the lavaan-formatted summary to console

## See also

[`SEMrun`](https://grassiMario.github.io/SEMgraph/reference/SEMrun.md).

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
sem1 <- SEMrun(sachs$graph, log(sachs$pkc), sachs$group, algo = "ricf")
#> RICF solver ended normally after 4 iterations 
#> 
#> deviance/df: 61.84434  srmr: 0.0700967 
#> 
#> Brown's combined P-value of node activation: 8.881784e-16 
#> 
#> Brown's combined P-value of node inhibition: 0.0004656197 
#> 
summary(sem1$fit)
#> RICF solver ended normally after 4 iterations 
#> 
#>   Estimator                                       ML 
#>   Optimization method                             RICF 
#> 
#>   Number of free parameters                       41 
#> 
#>   Number of observations                          1766 
#> 
#> Model Test User Model
#> 
#>   Test statistic (Deviance)                       2226.396 
#>   Degrees of freedom (df)                         36 
#>   Deviance/df                                     61.844 
#>   Standardized Root Mean Square Residual (srmr)   0.07 
#> 
#> Parameter Estimates:
#> 
#> Regressions: 
#> 
#>     lhs op   rhs    est
#> 1   Raf  ~ group -0.467
#> 2   Mek  ~ group  0.122
#> 3  Plcg  ~ group  0.268
#> 4  PIP2  ~ group  0.110
#> 5  PIP3  ~ group -0.301
#> 6   Erk  ~ group  0.179
#> 7   Akt  ~ group  0.150
#> 8   PKA  ~ group  0.290
#> 9   PKC  ~ group  0.116
#> 10  P38  ~ group -0.062
#> 11  Jnk  ~ group  0.333
#> 12  Mek  ~   Raf  0.736
#> 13  Erk  ~   Mek  0.016
#> 14 PIP2  ~  Plcg  0.199
#> 15  PKC  ~  Plcg -0.010
#> 16  PKC  ~  PIP2  0.029
#> 17 Plcg  ~  PIP3  0.288
#> 18 PIP2  ~  PIP3  0.528
#> 19  Akt  ~  PIP3 -0.023
#> 20  Raf  ~   PKA  0.008
#> 21  Mek  ~   PKA -0.017
#> 22  Erk  ~   PKA  0.407
#> 23  Akt  ~   PKA  0.507
#> 24  P38  ~   PKA  0.000
#> 25  Jnk  ~   PKA -0.002
#> 26  Raf  ~   PKC  0.007
#> 27  Mek  ~   PKC -0.015
#> 28  P38  ~   PKC  0.652
#> 29  Jnk  ~   PKC  0.096
#> 
#> Covariances: 
#> 
#>    lhs op  rhs    est
#> 36 PKA ~~ PIP3 -0.016
#> 
#> Variances: 
#> 
#>      lhs op   rhs   est
#> 30 group ~~ group 1.000
#> 31   Raf ~~   Raf 0.785
#> 32   Mek ~~   Mek 0.524
#> 33  Plcg ~~  Plcg 0.892
#> 34  PIP2 ~~  PIP2 0.653
#> 35  PIP3 ~~  PIP3 0.909
#> 37   Erk ~~   Erk 0.762
#> 38   Akt ~~   Akt 0.672
#> 39   PKA ~~   PKA 0.916
#> 40   PKC ~~   PKC 0.986
#> 41   P38 ~~   P38 0.580
#> 42   Jnk ~~   Jnk 0.873
#> 
```
