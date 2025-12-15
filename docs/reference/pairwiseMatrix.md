# Pairwise plotting of multivariate data

Display a pairwise scatter plot of two datasets for a random selection
of variables. If the second dataset is not given, the function displays
a histogram with normal curve superposition.

## Usage

``` r
pairwiseMatrix(x, y = NULL, size = nrow(x), r = 4, c = 4, ...)
```

## Arguments

- x:

  A matrix or data.frame (n x p) of continuous data.

- y:

  A matrix or data.frame (n x q) of continuous data.

- size:

  number of rows to be sampled (default `size = nrow(x)`).

- r:

  number of rows of the plot layout (default `r = 4`).

- c:

  number of columns of the plot layout (default `c = 4`).

- ...:

  Currently ignored.

## Value

No return value

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
adjdata <- SEMbap(sachs$graph, log(sachs$pkc))$data
#> DAG conversion : TRUE
#> Bow-free covariances search. Use method: cggm ...
#> Number of bow-free covariances / df : 4 / 37 
#> Max parent set(S) / Sparsity idx(s) : 5 / 18 
#> Number of clusters / number of nodes: 3 / 6 
#> 
rawdata <- log(sachs$pkc)
pairwiseMatrix(adjdata, rawdata, size = 1000)

```
