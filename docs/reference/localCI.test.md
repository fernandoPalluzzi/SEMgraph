# Conditional Independence (CI) local tests of an acyclic graph

P-values of one minimal testable implication (with the smallest possible
conditioning set) is returned per missing edge given an acyclic graph
(DAG or BAP) using the function
[`impliedConditionalIndependencies`](https://rdrr.io/pkg/dagitty/man/impliedConditionalIndependencies.html)
plus the function
[`localTests`](https://rdrr.io/pkg/dagitty/man/localTests.html) from
package `dagitty`. Without assuming any particular dependence structure,
the p-values of every CI test, in a DAG (BAP), is then combined using
the Bonferroni’s statistic in an overall test of the fitted model, B =
K\*min(p1,...,pK), as reviewed in Vovk & Wang (2020).

## Usage

``` r
localCI.test(graph, data, bap = FALSE, limit = 100, verbose = TRUE, ...)
```

## Arguments

- graph:

  A directed graph as an igraph object.

- data:

  A data matrix with subjects as rows and variables as columns.

- bap:

  If TRUE, the input graph is trasformend in a BAP, if FALSE (defult)
  the input graph is reduced in a DAG.

- limit:

  An integer value corresponding to the size of the extracted acyclic
  graph. Beyond this limit, switch to Shipley's C-test (Shipley 2000) is
  enabled to reduce the computational burden. By default, `limit = 100`.

- verbose:

  If TRUE, LocalCI results will be showed to screen (default = TRUE).

- ...:

  Currently ignored.

## Value

A list of three objects: (i) "dag": the DAG used to perform the localCI
test (ii) "msep": the list of all m-separation tests over missing edges
in the input graph and (iii) "mtest":the overall Bonferroni's P-value.

## References

Vovk V, Wang R (2020). Combining p-values via averaging. Biometrika
107(4): 791-808. \<https://doi.org/10.1093/biomet/asaa027\>

Shipley B (2000). A new inferential test for path models based on DAGs.
Structural Equation Modeling, 7(2): 206-218.
\<https://doi.org/10.1207/S15328007SEM0702_4\>

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
# Nonparanormal(npn) transformation
als.npn <- transformData(alsData$exprs)$data
#> Conducting the nonparanormal transformation via shrunkun ECDF...done.

sem <- SEMrun(alsData$graph, als.npn)
#> NLMINB solver ended normally after 1 iterations 
#> 
#> deviance/df: 10.92504  srmr: 0.2858859 
#> 
B_test <- localCI.test(sem$graph, als.npn, verbose = TRUE)
#> d-separation test (minimal set) of 420 edges...
#>     B_test df B_pvalue
#> 1 376.6275  1        0
```
