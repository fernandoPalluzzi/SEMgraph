# Parameter Estimates of a fitted SEM

Wrapper of the lavaan parameterEstimates() function for RICF and CGGM
algorithms

## Usage

``` r
parameterEstimates(fit, ...)
```

## Arguments

- fit:

  A RICF or constrained GGM fitted model object.

- ...:

  Currently ignored.

## Value

A data.frame containing the estimated parameters

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
ricf1 <- SEMrun(sachs$graph, log(sachs$pkc), sachs$group, algo = "ricf")
#> RICF solver ended normally after 4 iterations 
#> 
#> deviance/df: 61.84434  srmr: 0.0700967 
#> 
#> Brown's combined P-value of node activation: 8.881784e-16 
#> 
#> Brown's combined P-value of node inhibition: 0.0004656197 
#> 
parameterEstimates(ricf1$fit)
#>      lhs op   rhs    est
#> 1    Raf  ~ group -0.467
#> 2    Mek  ~ group  0.122
#> 3   Plcg  ~ group  0.268
#> 4   PIP2  ~ group  0.110
#> 5   PIP3  ~ group -0.301
#> 6    Erk  ~ group  0.179
#> 7    Akt  ~ group  0.150
#> 8    PKA  ~ group  0.290
#> 9    PKC  ~ group  0.116
#> 10   P38  ~ group -0.062
#> 11   Jnk  ~ group  0.333
#> 30 group ~~ group  1.000
#> 12   Mek  ~   Raf  0.736
#> 13   Erk  ~   Mek  0.016
#> 14  PIP2  ~  Plcg  0.199
#> 15   PKC  ~  Plcg -0.010
#> 16   PKC  ~  PIP2  0.029
#> 17  Plcg  ~  PIP3  0.288
#> 18  PIP2  ~  PIP3  0.528
#> 19   Akt  ~  PIP3 -0.023
#> 20   Raf  ~   PKA  0.008
#> 21   Mek  ~   PKA -0.017
#> 22   Erk  ~   PKA  0.407
#> 23   Akt  ~   PKA  0.507
#> 24   P38  ~   PKA  0.000
#> 25   Jnk  ~   PKA -0.002
#> 26   Raf  ~   PKC  0.007
#> 27   Mek  ~   PKC -0.015
#> 28   P38  ~   PKC  0.652
#> 29   Jnk  ~   PKC  0.096
#> 31   Raf ~~   Raf  0.785
#> 32   Mek ~~   Mek  0.524
#> 33  Plcg ~~  Plcg  0.892
#> 34  PIP2 ~~  PIP2  0.653
#> 35  PIP3 ~~  PIP3  0.909
#> 36   PKA ~~  PIP3 -0.016
#> 37   Erk ~~   Erk  0.762
#> 38   Akt ~~   Akt  0.672
#> 39   PKA ~~   PKA  0.916
#> 40   PKC ~~   PKC  0.986
#> 41   P38 ~~   P38  0.580
#> 42   Jnk ~~   Jnk  0.873

cggm1 <- SEMrun(sachs$graph, log(sachs$pkc), sachs$group, algo = "cggm")
#> DAG conversion : TRUE
#> GGM (de-biased nodewise L1) solver ended normally after 11 iterations 
#> 
#> deviance/df: 60.18925  srmr: 0.0703002 
#> 
#> Brown's combined P-value of node activation: 6.724732e-12 
#> 
#> Brown's combined P-value of node inhibition: 0.8709446 
#> 
parameterEstimates(cggm1$fit)
#>     lhs op   rhs    est    se  z_test pvalue ci.lower ci.uppper
#> 1   PKA  ~ group  0.290 0.023  12.720  0.000    0.245     0.334
#> 2   PKC  ~ group  0.115 0.024   4.766  0.000    0.068     0.162
#> 5  PIP3  ~ group -0.301 0.023 -13.278  0.000   -0.346    -0.257
#> 6   Mek  ~ group  0.120 0.020   5.929  0.000    0.080     0.159
#> 10  Raf  ~ group -0.465 0.022 -21.032  0.000   -0.509    -0.422
#> 13 PIP2  ~ group  0.108 0.021   5.188  0.000    0.067     0.149
#> 16 Plcg  ~ group  0.266 0.024  11.308  0.000    0.220     0.312
#> 18  Jnk  ~ group  0.332 0.023  14.237  0.000    0.286     0.378
#> 21  P38  ~ group -0.061 0.019  -3.223  0.001   -0.099    -0.024
#> 24  Akt  ~ group  0.149 0.021   7.040  0.000    0.108     0.191
#> 27  Erk  ~ group  0.178 0.022   8.007  0.000    0.134     0.221
#> 3   PKC  ~  PIP2  0.028 0.025   1.101  0.271   -0.022     0.077
#> 4   PKC  ~  Plcg -0.008 0.025  -0.314  0.753   -0.058     0.042
#> 7   Mek  ~   PKA -0.016 0.018  -0.884  0.376   -0.051     0.019
#> 8   Mek  ~   PKC -0.014 0.017  -0.823  0.410   -0.048     0.020
#> 9   Mek  ~   Raf  0.734 0.019  37.816  0.000    0.696     0.772
#> 11  Raf  ~   PKA  0.007 0.022   0.303  0.762   -0.036     0.050
#> 12  Raf  ~   PKC  0.006 0.021   0.297  0.767   -0.035     0.048
#> 14 PIP2  ~  PIP3  0.527 0.021  25.119  0.000    0.486     0.568
#> 15 PIP2  ~  Plcg  0.198 0.020   9.751  0.000    0.158     0.238
#> 17 Plcg  ~  PIP3  0.286 0.024  12.159  0.000    0.240     0.332
#> 19  Jnk  ~   PKA  0.000 0.023  -0.018  0.985   -0.046     0.045
#> 20  Jnk  ~   PKC  0.095 0.022   4.245  0.000    0.051     0.139
#> 22  P38  ~   PKA  0.000 0.019   0.000  1.000   -0.037     0.037
#> 23  P38  ~   PKC  0.651 0.018  35.721  0.000    0.616     0.687
#> 25  Akt  ~   PKA  0.506 0.020  24.844  0.000    0.466     0.546
#> 26  Akt  ~  PIP3 -0.023 0.020  -1.112  0.266   -0.063     0.017
#> 28  Erk  ~   PKA  0.406 0.022  18.741  0.000    0.364     0.449
#> 29  Erk  ~   Mek  0.015 0.021   0.688  0.491   -0.027     0.056
```
