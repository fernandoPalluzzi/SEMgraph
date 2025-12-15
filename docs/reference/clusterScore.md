# Module scoring

Generate factor scores, principal component scores, or projection scores
of latent, composite, and unmeasured variable modules, respectively, and
fit them with an exogenous group effect.

## Usage

``` r
clusterScore(
  graph,
  data,
  group,
  HM = "LV",
  type = "wtc",
  size = 5,
  verbose = FALSE,
  ...
)
```

## Arguments

- graph:

  An igraph object.

- data:

  A matrix or data.frame. Rows correspond to subjects, and columns to
  graph nodes.

- group:

  A binary vector. This vector must be as long as the number of
  subjects. Each vector element must be 1 for cases and 0 for control
  subjects.

- HM:

  Hidden model type. For each defined hidden module: (i) if `HM = "LV"`,
  a latent variable (LV) will be defined as common unknown cause acting
  on cluster nodes; (ii) if `HM = "CV"`, cluster nodes will be
  considered as regressors of a latent composite variable (CV); (iii) if
  `HM = "UV"`, an unmeasured variable (UV) model will be generated for
  each module, where source nodes (i.e., in-degree = 0) act as common
  regressors influencing the other nodes via an unmeasured variable. By
  default, HM is set to "LV" (i.e., the latent variable model).

- type:

  Graph clustering method. If `type = "tahc"`, network modules are
  generated using the tree agglomerative hierarchical clustering method
  (Yu et al., 2015). Other non-tree clustering methods from igraph
  package include: "wtc" (default value; walktrap community structure
  with short random walks), "ebc" (edge betweenness clustering), "fgc"
  (fast greedy method), "lbc" (label propagation method), "lec" (leading
  eigenvector method), "loc" (multi-level optimization), "opc" (optimal
  communiy structure), "sgc" (spinglass statistical mechanics). By
  default, the "wtc" method is used.

- size:

  Minimum number of nodes per hidden module. By default, a minimum
  number of 5 nodes is required.

- verbose:

  A logical value. If TRUE, intermediate graphs will be displayed during
  the execution. In addition, a reduced graph with clusters as nodes
  will be fitted and showed to screen (see also
  [`mergeNodes`](https://grassiMario.github.io/SEMgraph/reference/mergeNodes.md)).
  By default, `verbode = FALSE`.

- ...:

  Currently ignored.

## Value

A list of 3 objects:

1.  "fit", hidden module fitting as a lavaan object;

2.  "membership", hidden module nodes membership;
    [`clusterGraph`](https://grassiMario.github.io/SEMgraph/reference/clusterGraph.md)
    function;

3.  "dataHM", data matrix with cluster scores in first columns.

## References

Grassi M, Palluzzi F, Tarantino B (2022). SEMgraph: An R Package for
Causal Network Analysis of High-Throughput Data with Structural Equation
Models. Bioinformatics, 38 (20), 4829–4830
\<https://doi.org/10.1093/bioinformatics/btac567\>

## See also

See
[`clusterGraph`](https://grassiMario.github.io/SEMgraph/reference/clusterGraph.md)
and [`cplot`](https://grassiMario.github.io/SEMgraph/reference/cplot.md)
for graph clustering.

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
# Nonparanormal(npn) transformation
als.npn <- transformData(alsData$exprs)$data
#> Conducting the nonparanormal transformation via shrunkun ECDF...done.

C <- clusterScore(graph = alsData$graph, data = als.npn,
                  group = alsData$group,
                  HM = "LV",
                  type = "wtc",
                  verbose = FALSE)
#> modularity = 0.528642 
#> 
#> Community sizes
#>  4  1  3  2 
#>  3  8  9 11 
#> 
#> Model converged: TRUE 
#> SRMR: 6.89267e-10 
#> 
summary(C$fit)
#> lavaan 0.6-20 ended normally after 16 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         9
#> 
#>   Number of observations                           160
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 0
#> 
#> Parameter Estimates:
#> 
#>   Standard errors                             Standard
#>   Information                                 Expected
#>   Information saturated (h1) model          Structured
#> 
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   LV1 ~                                               
#>     group            -0.848    0.251   -3.377    0.001
#>   LV2 ~                                               
#>     group             0.869    0.226    3.847    0.000
#>   LV3 ~                                               
#>     group            -0.987    0.225   -4.382    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>  .LV1 ~~                                              
#>    .LV2              -0.201    0.083   -2.416    0.016
#>    .LV3               0.787    0.103    7.667    0.000
#>  .LV2 ~~                                              
#>    .LV3              -0.367    0.079   -4.650    0.000
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .LV1               1.152    0.129    8.944    0.000
#>    .LV2               0.930    0.104    8.944    0.000
#>    .LV3               0.925    0.103    8.944    0.000
#> 
head(C$dataHM)
#>      group        LV1        LV2         LV3       317        572        581
#> ALS2     1 -0.7989464 -0.7527118 -0.42868026 1.0661487 -0.7577341 -1.1225852
#> ALS3     1 -0.8635143 -0.8612364 -0.68134940 1.6294444 -2.1342965 -0.2662029
#> ALS4     1  0.9117982 -1.1961490  0.73949684 0.9870911  0.3160400 -0.2006991
#> ALS5     1 -1.1279817 -0.2024923 -1.35460499 0.5799344 -0.9135118 -1.0127737
#> ALS6     1 -0.2946798  1.3716185 -0.07552589 0.1360061 -0.8900137 -1.4346740
#> ALS7     1 -1.7251688  0.4125732 -1.46126288 0.5987578 -1.6885942 -0.4181436
#>             596          598         836        842      54205        1616
#> ALS2  0.3666365 -0.007977011 -0.08785375 -0.8003346  0.2497352 -0.20069907
#> ALS3  1.3168559 -0.716404714  0.89001373 -0.1199272  0.6370285 -0.77886604
#> ALS4  1.2812248 -0.216987905 -0.48853105  0.9375004 -1.1521227  0.24973518
#> ALS5  1.6885942 -0.418143644  1.21426434 -0.2827396  0.8443656 -0.15211855
#> ALS6  0.6370285  0.656498936  0.10387774 -0.1038777 -0.2497352  0.05586618
#> ALS7 -0.6177855 -2.299257173  0.63702853 -1.2470790  1.2812248 -2.29925717
#>           79139        5606       5608        4217       5600       6300
#> ALS2  0.9870911 -0.03989473  0.8003346 -0.18446080 -1.0661487 -1.2142643
#> ALS3  1.3168559 -1.62944439  1.3933390 -0.03989473 -0.8900137 -0.8900137
#> ALS4 -0.8003346  0.28273958 -0.3666365  0.69617310  0.2827396  0.9870911
#> ALS5  2.0111279  0.02393297  0.9620141  1.57527382 -0.6177855 -2.5616910
#> ALS6 -0.2169879  0.10387774  0.5613048  0.75773406 -1.2470790 -1.3168559
#> ALS7  1.0127737 -0.57993445  2.2992572  1.21426434  0.1199272  0.4008634
#>              5603      1432        4744       4747       4741       5530
#> ALS2  0.332813985 0.8669746 -0.67620930 -0.5799344 -0.7577341  0.6961731
#> ALS3 -0.103877738 0.2497352 -0.89001373 -0.8003346 -0.8443656  0.8443656
#> ALS4  0.266202856 1.1225852 -1.47848512 -1.5251790 -1.1521227 -0.5064807
#> ALS5  1.316855924 2.1342965  0.15211855  0.2333318 -0.2497352  0.9870911
#> ALS6 -0.233331753 1.3933390  1.57527382  0.7577341  1.3541545 -0.8669746
#> ALS7 -0.007977011 1.0127737 -0.05586618  0.5245871  0.4181436  2.5616910
#>             5532       5533       5534        5535        5630       6647
#> ALS2  0.41814364  0.4530701  0.2993503  0.03989473 -1.31685592 -0.7369195
#> ALS3  0.57993445 -1.2142643  0.5987578  1.03910878 -0.54285880  0.7577341
#> ALS4 -0.73691953 -0.3666365 -0.6564989  1.62944439  1.18264983 -1.7539739
#> ALS5  1.21426434  0.3328140  1.4346740  0.71640471 -0.98709115  0.3328140
#> ALS6  0.07185121 -0.1199272  0.3836965 -0.47073008 -0.40086343  1.2812248
#> ALS7  1.06614870 -0.6177855  1.3168559  0.43554360 -0.03989473 -0.1844608
#>             7132       7133      10452      84134
#> ALS2 -0.07185121 -0.7788660 -0.9620141 -0.1199272
#> ALS3 -1.12258519 -0.9620141 -1.1521227 -1.1826498
#> ALS4  0.40086343  0.1038777 -1.2142643 -0.2993503
#> ALS5  0.20069907 -2.2992572 -1.9115966 -0.8669746
#> ALS6  0.91351183 -0.6961731 -1.6885942 -0.5245871
#> ALS7  0.96201413  0.8003346 -0.5799344 -0.1682687
table(C$membership)
#> 
#>  1  2  3 
#>  8 11  9 
```
