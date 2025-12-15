# Graph weighting methods

Add data-driven edge and node weights to the input graph.

## Usage

``` r
weightGraph(graph, data, group = NULL, method = "r2z", limit = 10000, ...)
```

## Arguments

- graph:

  An igraph object.

- data:

  A matrix or data.frame. Rows correspond to subjects, and columns to
  graph nodes.

- group:

  Binary vector. This vector must be as long as the number of subjects.
  Each vector element must be 1 for cases and 0 for control subjects. By
  default, `group = NULL`. If group is not NULL, also node weighting is
  actived, and node weights correspond to the attribute: V(graph)\$pv
  (P-value of the z-test = b/SE(b) from simple linear regression y ~ x,
  i.e., lm(node ~ group)) and V(graph)\$sign (-1 if z\<-2, +1 if z\>2, 0
  otherwise).

- method:

  Edge weighting method. It can be one of the following:

  1.  "r2z", weight edges are defined using Fisher's r-to-z transform
      (Fisher, 1915) to test the correlation coefficient of pairs of
      interacting nodes, if `group=NULL`. Otherwise, the difference
      between group of the r-to-z trasform will be tested. Edge weights
      correspond to the attribute: E(graph)\$pv (P-value of the z-test)
      and E(graph)\$sign (-1 if z\<-2, +1 if z\>2, 0 otherwise).

  2.  "sem", edge weights are defined by a SEM model that implies
      testing the group effect simultaneously on source and sink nodes.
      A new parameter w is defined as the weighted sum of the total
      effect of the group on source and sink nodes, adjusted by node
      degree centrality. Edge weights correspond to the attribute:
      E(graph)\$pv (P-value of the z-test = w/SE(w)) and E(graph)\$sign
      (-1 if z\<-2, +1 if z\>2, 0 otherwise). Not available if
      `group=NULL`.

  3.  "cov", edge weights are defined by a new parameter w combining the
      group effect on the source node (mean group difference, adjusted
      by source degree centrality), the sink node (mean group
      difference, adjusted by sink degree centrality), and the
      source–sink interaction (correlation difference). Edge weights
      correspond to the attribute: E(graph)\$pv (P-value of the z-test =
      w/SE(w) of the combined difference of the group over source node,
      sink node, and their connection) and E(graph)\$sign (-1 if z\<-2,
      +1 if z\>2, 0 otherwise). Not available if `group=NULL`.

  4.  "cfa", edge weights are defined by a CFA1 model that implies
      testing the group effect, w on a latent variable (LV) with
      observed indicators two interacting nodes, fixing loading
      coefficients and residual variances for model identification. Edge
      weights correspond to the attribute: E(graph)\$pv (P-value of the
      z-test = w/SE(w) of the group effect on the LV) and E(graph)\$sign
      (-1 if z\<-2, +1 if z\>2, 0 otherwise). Not available if
      `group=NULL`.

- limit:

  An integer value corresponding to the number of graph edges. Beyond
  this limit, multicore computation is enabled to reduce the
  computational burden. By default, `limit = 10000`.

- ...:

  Currently ignored.

## Value

A weighted graph, as an igraph object.

## References

Grassi M, Tarantino B (2023). \[Supplementary material of\] SEMtree:
tree-based structure learning methods with structural equation models.
Bioinformatics, 39 (6), 4829–4830
\<https://doi.org/10.1093/bioinformatics/btad377\>

Fisher RA (1915). Frequency Distribution of the Values of the
Correlation Coefficient in Samples from an Indefinitely Large
Population. Biometrika, 10(4), 507–521. \<doi:10.2307/2331838\>

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
# Graph weighting
G <- weightGraph(graph = sachs$graph,
                 data = log(sachs$pkc),
                 group = sachs$group,
                 method = "r2z")

# New edge attributes
head(E(G)$pv); summary(E(G)$pv)
#> [1] 2.661577e-01 8.370494e-01 5.528435e-01 4.401098e-01 4.742076e-01
#> [6] 1.003647e-09
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.000e+00 4.400e-07 2.413e-01 2.953e-01 5.031e-01 9.643e-01 
head(E(G)$zsign); table(E(G)$zsign)
#> [1] 0 0 0 0 0 1
#> 
#>  0  1 
#> 13  7 

# New node attributes
head(V(G)$pv); summary(V(G)$pv)
#> [1] 1.721660e-35 1.576567e-06 2.219336e-38 7.635299e-22 8.952626e-95
#> [6] 5.688104e-01
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.000e+00 0.000e+00 0.000e+00 1.076e-01 7.900e-07 6.152e-01 
head(V(G)$zsign); table(V(G)$zsign)
#> [1]  1  1 -1 -1 -1  0
#> 
#> -1  0  1 
#>  3  2  6 
```
