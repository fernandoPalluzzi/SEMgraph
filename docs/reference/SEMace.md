# Compute the Average Causal Effect (ACE) for a given source-sink pair

Compute total effects as ACEs of variables X on variables Y in a
directed acyclic graph (DAG). The ACE will be estimated as the path
coefficient of X (i.e., theta) in the linear equation Y ~ X + Z. The set
Z is defined as the adjustment (or conditioning) set of Y over X,
applying various adjustement sets. Standard errors (SE), for each ACE,
are computed following the `lm` standard procedure or a bootstrap-based
procedure (see [`boot`](https://rdrr.io/pkg/boot/man/boot.html) for
details).

## Usage

``` r
SEMace(
  graph,
  data,
  group = NULL,
  type = "parents",
  effect = "all",
  method = "BH",
  alpha = 0.05,
  boot = NULL,
  ...
)
```

## Arguments

- graph:

  An igraph object.

- data:

  A matrix or data.frame. Rows correspond to subjects, and columns to
  graph nodes (variables).

- group:

  A binary vector. This vector must be as long as the number of
  subjects. Each vector element must be 1 for cases and 0 for control
  subjects. If `group = NULL` (default), group influence will not be
  considered.

- type:

  character Conditioning set Z. If "parents" (default) the Pearl's
  back-door set (Pearl, 1998), "minimal" the dagitty minimal set
  (Perkovic et al, 2018), or "optimal" the O-set with the smallest
  asymptotic variance (Witte et al, 2020) are computed.

- effect:

  character X to Y effect. If "all" (default) all effects from X to Y,
  "source2sink" only effects from source X to sink Y, or "direct" only
  direct effects from X to Y are computed.

- method:

  Multiple testing correction method. One of the values available in
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). By default,
  `method = "BH"` (i.e., FDR multiple test correction).

- alpha:

  Significance level for ACE selection (by default, `alpha = 0.05`).

- boot:

  The number of bootstrap samplings enabling bootstrap computation of
  ACE standard errors. If `NULL` (default), bootstrap is disabled.

- ...:

  Currently ignored.

## Value

A data.frame of ACE estimates between network sources and sinks.

## References

Pearl J (1998). Graphs, Causality, and Structural Equation Models.
Sociological Methods & Research, 27(2):226-284.
\<https://doi.org/10.1177/0049124198027002004\>

Perkovic E, Textor J, Kalisch M, Maathuis MH (2018). Complete graphical
characterization and construction of adjustment sets in Markov
equivalence classes of ancestral graphs. Journal of Machine Learning
Research, 18:1-62. \<http://jmlr.org/papers/v18/16-319.html\>

Witte J, Henckel L, Maathuis MH, Didelez V (2020). On efficient
adjustment in causal graphs. Journal of Machine Learning Research,
21:1-45. \<http://jmlr.org/papers/v21/20-175.htm\>

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
# ACE without group, O-set, all effects:
ace1 <- SEMace(graph = sachs$graph, data = log(sachs$pkc),
               group = NULL, type = "optimal", effect = "all",
               method = "BH", alpha = 0.05, boot = NULL)
#> 
#> WARNING: input graph is not acyclic!
#>  Applying graph -> DAG conversion.
#> DAG conversion : TRUE
#> 
#>  Frequency distribution of path length from X to Y :
#>  1  2  3  4 
#> 18 11  6  1 
#> 
#>  ACE= 1 of 36 ACE= 2 of 36 ACE= 3 of 36 ACE= 4 of 36 ACE= 5 of 36 ACE= 6 of 36 ACE= 7 of 36 ACE= 8 of 36 ACE= 9 of 36 ACE= 10 of 36 ACE= 11 of 36 ACE= 12 of 36 ACE= 13 of 36 ACE= 14 of 36 ACE= 15 of 36 ACE= 16 of 36 ACE= 17 of 36 ACE= 18 of 36 ACE= 19 of 36 ACE= 20 of 36 ACE= 21 of 36 ACE= 22 of 36 ACE= 23 of 36 ACE= 24 of 36 ACE= 25 of 36 ACE= 26 of 36 ACE= 27 of 36 ACE= 28 of 36 ACE= 29 of 36 ACE= 30 of 36 ACE= 31 of 36 ACE= 32 of 36 ACE= 33 of 36 ACE= 34 of 36 ACE= 35 of 36 ACE= 36 of 36
print(ace1)
#>       pathL sink op source    est    se      z pvalue ci.lower ci.upper
#> PKA       1  Mek <-    PKA -0.074 0.024 -3.119  0.002   -0.121   -0.028
#> PKA1      1  Raf <-    PKA -0.124 0.024 -5.258  0.000   -0.171   -0.078
#> PKA2      1  Jnk <-    PKA  0.093 0.024  3.938  0.000    0.047    0.139
#> PKA4      1  Akt <-    PKA  0.546 0.020 27.439  0.000    0.507    0.585
#> PKA5      1  Erk <-    PKA  0.458 0.021 21.582  0.000    0.416    0.499
#> PKC2      1  Jnk <-    PKC  0.128 0.024  5.436  0.000    0.082    0.174
#> PKC3      1  P38 <-    PKC  0.646 0.018 35.484  0.000    0.611    0.682
#> PIP3      3  Mek <-   PIP3  0.085 0.024  3.560  0.000    0.038    0.131
#> PIP31     3  Raf <-   PIP3  0.155 0.023  6.600  0.000    0.109    0.201
#> 11        1 PIP2 <-   PIP3  0.536 0.020 26.687  0.000    0.497    0.576
#> 12        1 Plcg <-   PIP3  0.207 0.023  8.892  0.000    0.161    0.253
#> PIP32     3  Jnk <-   PIP3 -0.077 0.024 -3.250  0.001   -0.124   -0.031
#> PIP34     1  Akt <-   PIP3 -0.065 0.020 -3.249  0.001   -0.104   -0.026
#> PIP35     4  Erk <-   PIP3 -0.063 0.021 -2.988  0.003   -0.105   -0.022
#> Raf       1  Mek <-    Raf  0.683 0.018 38.898  0.000    0.649    0.718
#> Raf1      2  Erk <-    Raf -0.081 0.021 -3.784  0.000   -0.122   -0.039
#> Plcg1     2  Mek <-   Plcg -0.059 0.024 -2.451  0.014   -0.107   -0.012
#> Plcg2     2  Raf <-   Plcg -0.091 0.024 -3.822  0.000   -0.138   -0.044
#> Plcg3     1 PIP2 <-   Plcg  0.227 0.020 11.435  0.000    0.188    0.265
#> Plcg4     2  Jnk <-   Plcg  0.104 0.024  4.317  0.000    0.057    0.151
#> Plcg6     3  Erk <-   Plcg  0.057 0.022  2.631  0.009    0.015    0.099

# ACE with group perturbation, Pa-set, direct effects:
ace2 <- SEMace(graph = sachs$graph, data = log(sachs$pkc),
               group = sachs$group, type = "parents", effect = "direct",
               method = "none", alpha = 0.05, boot = NULL)
#> 
#>  Frequency distribution of path length from X to Y :
#>  1 
#> 20 
#> 
#>  ACE= 1 of 20 ACE= 2 of 20 ACE= 3 of 20 ACE= 4 of 20 ACE= 5 of 20 ACE= 6 of 20 ACE= 7 of 20 ACE= 8 of 20 ACE= 9 of 20 ACE= 10 of 20 ACE= 11 of 20 ACE= 12 of 20 ACE= 13 of 20 ACE= 14 of 20 ACE= 15 of 20 ACE= 16 of 20 ACE= 17 of 20 ACE= 18 of 20 ACE= 19 of 20 ACE= 20 of 20
print(ace2)
#>       pathL sink op source d_est  d_se    d_z pvalue d_lower d_upper
#> PKA       1 PIP3 <-    PKA 0.000 0.000  3.584  0.000   0.000   0.000
#> PKA5      1  Akt <-    PKA 0.341 0.039  8.722  0.000   0.265   0.418
#> PKA6      1  Erk <-    PKA 0.242 0.042  5.706  0.000   0.159   0.326
#> PKC2      1  Jnk <-    PKC 0.643 0.044 14.506  0.000   0.556   0.730
#> PKC3      1  P38 <-    PKC 0.345 0.036  9.684  0.000   0.275   0.415
#> PIP31     1 PIP2 <-   PIP3 0.351 0.040  8.693  0.000   0.272   0.431
#> PIP32     1 Plcg <-   PIP3 0.321 0.047  6.774  0.000   0.228   0.414
#> PIP33     1  Akt <-   PIP3 0.079 0.040  1.994  0.046   0.001   0.158
#> Plcg1     1 PIP2 <-   Plcg 0.375 0.039  9.715  0.000   0.299   0.450
```
