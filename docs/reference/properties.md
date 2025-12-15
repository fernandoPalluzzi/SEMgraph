# Graph properties summary and graph decomposition

Produces a summary of network properties and returns graph components
(ordered by decreasing size), without self-loops.

## Usage

``` r
properties(graph, data = NULL, ...)
```

## Arguments

- graph:

  Input network as an igraph object.

- data:

  An optional data matrix (default data = NULL) whith rows corresponding
  to subjects, and columns to graph nodes (variables). Nodes will be
  mapped onto variable names.

- ...:

  Currently ignored.

## Value

List of graph components, ordered by decreasing size (the first
component is the giant one), without self-loops.

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
# Extract the "Neurotrophin signaling pathway":
g <- kegg.pathways[["Neurotrophin signaling pathway"]]
summary(g)
#> IGRAPH 1033343 DNW- 119 357 -- 
#> + attr: name (v/c), weight (e/n)
properties(g)
#> Frequency distribution of graph components
#> 
#>   n.nodes n.graphs
#> 1     117        1
#> 
#> Percent of vertices in the giant component: 98.3 %
#> 
#>   is.simple      is.dag is.directed is.weighted 
#>        TRUE       FALSE        TRUE        TRUE 
#> 
#> which.mutual.FALSE 
#>                357 
#> [[1]]
#> IGRAPH e8455ea DNW- 117 357 -- 
#> + attr: name (v/c), weight (e/n)
#> + edges from e8455ea (vertex names):
#>  [1] 4914->5335   4914->5336   4914->25759  4914->399694 4914->53358 
#>  [6] 4914->6464   4914->2885   4914->10818  4914->57498  4914->4145  
#> [11] 4914->3667   4914->10603  4914->10019  4914->25970  4914->25    
#> [16] 4915->5335   4915->5336   4915->25759  4915->399694 4915->53358 
#> [21] 4915->6464   4915->2885   4915->4145   4915->3667   4915->10603 
#> [26] 4915->10019  4915->25970  4915->25     4916->5335   4916->5336  
#> [31] 4916->25759  4916->399694 4916->53358  4916->6464   4916->2885  
#> [36] 4916->4145   4916->3667   4916->10603  4916->10019  4916->25970 
#> + ... omitted several edges
#> 
```
