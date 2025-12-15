# GGM model summary

Generate a summary for a constrained Gaussian Graphical Model (GGM)
similar to lavaan-formated summary

## Usage

``` r
# S3 method for class 'GGM'
summary(object, ...)
```

## Arguments

- object:

  A constrained GGM fitted model object.

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
sem0 <- SEMrun(sachs$graph, log(sachs$pkc), algo = "cggm")
#> DAG conversion : TRUE
#> GGM (de-biased nodewise L1) solver ended normally after 9 iterations 
#> 
#> deviance/df: 66.94042  srmr: 0.0896276 
#> 
summary(sem0$fit)
#> GGM (de-biased nodewise L1) solver ended normally after 0 iterations 
#> 
#>   Estimator                                       ML 
#>   Optimization method                             GGM 
#> 
#>   Number of free parameters                       27 
#> 
#>   Number of observations                          1766 
#> 
#> Model Test User Model
#> 
#>   Test statistic (Deviance)                       2409.855 
#>   Degrees of freedom (df)                         36 
#>   Deviance/df                                     66.94 
#>   Standardized Root Mean Square Residual (srmr)   0.09 
#> 
#> Parameter Estimates:
#> 
#> Regressions:
#> 
#>     lhs op  rhs    est    se z_test pvalue ci.lower ci.uppper
#> 1   PKC  ~ PIP2  0.019 0.025  0.753  0.451   -0.030     0.068
#> 2   PKC  ~ Plcg  0.014 0.025  0.541  0.588   -0.036     0.063
#> 3   Mek  ~  PKA  0.010 0.018  0.547  0.584   -0.025     0.044
#> 4   Mek  ~  PKC -0.005 0.017 -0.262  0.793   -0.039     0.030
#> 5   Mek  ~  Raf  0.682 0.018 38.894  0.000    0.648     0.717
#> 6   Raf  ~  PKA -0.123 0.024 -5.225  0.000   -0.170    -0.077
#> 7   Raf  ~  PKC -0.037 0.024 -1.548  0.122   -0.083     0.010
#> 8  PIP2  ~ PIP3  0.488 0.020 24.675  0.000    0.450     0.527
#> 9  PIP2  ~ Plcg  0.226 0.020 11.408  0.000    0.187     0.265
#> 10 Plcg  ~ PIP3  0.207 0.023  8.895  0.000    0.161     0.253
#> 11  Jnk  ~  PKA  0.092 0.024  3.903  0.000    0.046     0.138
#> 12  Jnk  ~  PKC  0.127 0.024  5.402  0.000    0.081     0.173
#> 13  P38  ~  PKA -0.017 0.018 -0.930  0.352   -0.053     0.019
#> 14  P38  ~  PKC  0.645 0.018 35.465  0.000    0.610     0.681
#> 15  Akt  ~  PKA  0.545 0.020 27.426  0.000    0.506     0.584
#> 16  Akt  ~ PIP3 -0.064 0.020 -3.208  0.001   -0.103    -0.025
#> 17  Erk  ~  PKA  0.455 0.021 21.475  0.000    0.414     0.497
#> 18  Erk  ~  Mek -0.020 0.021 -0.931  0.352   -0.061     0.022
#> 
#> Variances:
#> 
#>   PKA   PKC  PIP3   Mek   Raf  PIP2  Plcg   Jnk   P38   Akt   Erk 
#> 1.000 1.000 1.000 0.535 0.983 0.663 0.957 0.974 0.583 0.691 0.790 
```
