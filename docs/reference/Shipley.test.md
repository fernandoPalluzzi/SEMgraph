# Missing edge testing implied by a DAG with Shipley's basis-set

Compute all the P-values of the d-separation tests implied by the
missing edges of a given acyclic graph (DAG). The conditioning set Z is
represented, in a DAG, by the union of the parent sets of X and Y
(Shipley, 2000). The results of every test, in a DAG, is then combined
using the Fisher’s statistic in an overall test of the fitted model C =
-2\*sum(log(P-value(k))), where C is distributed as a chi-squared
variate with df = 2k, as suggested by Shipley (2000).

## Usage

``` r
Shipley.test(
  graph,
  data,
  MCX2 = FALSE,
  cmax = Inf,
  limit = 100,
  verbose = TRUE,
  ...
)
```

## Arguments

- graph:

  A directed graph as an igraph object.

- data:

  A data matrix with subjects as rows and variables as columns.

- MCX2:

  If TRUE, a Monte Carlo P-value of the combined C test is enabled using
  the R code of Shipley extracted from
  \<https://github.com/BillShipley/CauseAndCorrelation\>.

- cmax:

  Maximum number of parents set, C. This parameter can be used to
  perform only those tests where the number of conditioning variables
  does not exceed the given value. High-dimensional conditional
  independence tests can be very unreliable. By default, cmax = Inf.

- limit:

  An integer value corresponding to the graph size (vcount) tolerance.
  Beyond this limit, multicore computation is enabled to reduce the
  computational burden. By default, `limit = 100`.

- verbose:

  If TRUE, Shipley's test results will be showed to screen (default =
  TRUE).

- ...:

  Currently ignored.

## Value

A list of three objects: (i) "dag": the DAG used to perform the Shipley
test (ii) "dsep": the data.frame of all d-separation tests over missing
edges in the DAG and (iii) "ctest": the overall Shipley's' P-value.

## References

Shipley B (2000). A new inferential test for path models based on DAGs.
Structural Equation Modeling, 7(2): 206-218.
\<https://doi.org/10.1207/S15328007SEM0702_4\>

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
#\donttest{
# Nonparanormal(npn) transformation
als.npn <- transformData(alsData$exprs)$data
#> Conducting the nonparanormal transformation via shrunkun ECDF...done.

sem <- SEMrun(alsData$graph, als.npn)
#> NLMINB solver ended normally after 1 iterations 
#> 
#> deviance/df: 10.92504  srmr: 0.2858859 
#> 
C_test <- Shipley.test(sem$graph, als.npn, MCX2 = FALSE)
#> d-separation test (basis set) of 420 edges...
#>  / 0 % elapsed=00s    - 1 % elapsed=00s, remaining~00s \ 2 % elapsed=00s, remaining~00s | 4 % elapsed=00s, remaining~00s / 5 % elapsed=00s, remaining~00s - 6 % elapsed=00s, remaining~00s \ 7 % elapsed=00s, remaining~00s | 8 % elapsed=00s, remaining~00s / 10% elapsed=00s, remaining~00s - 11% elapsed=00s, remaining~00s \ 12% elapsed=00s, remaining~00s | 13% elapsed=00s, remaining~00s / 14% elapsed=00s, remaining~00s - 15% elapsed=00s, remaining~00s \ 17% elapsed=00s, remaining~00s | 18% elapsed=00s, remaining~00s / 19% elapsed=00s, remaining~00s - 20% elapsed=00s, remaining~00s \ 21% elapsed=00s, remaining~00s | 23% elapsed=00s, remaining~00s / 24% elapsed=00s, remaining~00s - 25% elapsed=00s, remaining~00s \ 26% elapsed=00s, remaining~00s | 27% elapsed=00s, remaining~00s / 29% elapsed=00s, remaining~00s - 30% elapsed=00s, remaining~00s \ 31% elapsed=00s, remaining~00s | 32% elapsed=00s, remaining~00s / 33% elapsed=00s, remaining~00s - 35% elapsed=00s, remaining~00s \ 36% elapsed=00s, remaining~00s | 37% elapsed=00s, remaining~00s / 38% elapsed=00s, remaining~00s - 39% elapsed=00s, remaining~00s \ 40% elapsed=00s, remaining~00s | 42% elapsed=00s, remaining~00s / 43% elapsed=00s, remaining~00s - 44% elapsed=00s, remaining~00s \ 45% elapsed=00s, remaining~00s | 46% elapsed=00s, remaining~00s / 48% elapsed=00s, remaining~00s - 49% elapsed=00s, remaining~00s \ 50% elapsed=00s, remaining~00s | 51% elapsed=00s, remaining~00s / 52% elapsed=00s, remaining~00s - 54% elapsed=00s, remaining~00s \ 55% elapsed=00s, remaining~00s | 56% elapsed=00s, remaining~00s / 57% elapsed=00s, remaining~00s - 58% elapsed=00s, remaining~00s \ 60% elapsed=00s, remaining~00s | 61% elapsed=00s, remaining~00s / 62% elapsed=00s, remaining~00s - 63% elapsed=00s, remaining~00s \ 64% elapsed=00s, remaining~00s | 65% elapsed=00s, remaining~00s / 67% elapsed=00s, remaining~00s - 68% elapsed=00s, remaining~00s \ 69% elapsed=00s, remaining~00s | 70% elapsed=00s, remaining~00s / 71% elapsed=00s, remaining~00s - 73% elapsed=00s, remaining~00s \ 74% elapsed=00s, remaining~00s | 75% elapsed=00s, remaining~00s / 76% elapsed=00s, remaining~00s - 77% elapsed=00s, remaining~00s \ 79% elapsed=00s, remaining~00s | 80% elapsed=00s, remaining~00s / 81% elapsed=00s, remaining~00s - 82% elapsed=00s, remaining~00s \ 83% elapsed=00s, remaining~00s | 85% elapsed=00s, remaining~00s / 86% elapsed=00s, remaining~00s - 87% elapsed=00s, remaining~00s \ 88% elapsed=00s, remaining~00s | 89% elapsed=00s, remaining~00s / 90% elapsed=00s, remaining~00s - 92% elapsed=00s, remaining~00s \ 93% elapsed=00s, remaining~00s | 94% elapsed=00s, remaining~00s / 95% elapsed=00s, remaining~00s - 96% elapsed=00s, remaining~00s \ 98% elapsed=00s, remaining~00s | 99% elapsed=00s, remaining~00s / 100% elapsed=00s, remaining~00s
#>     C_test  df pvalue
#> 1 7038.822 840      0
#MC_test <- Shipley.test(sem$graph, als.npn, MCX2 = TRUE)
#}
```
