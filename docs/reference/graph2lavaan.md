# Graph to lavaan model

Convert an igraph object to a model (lavaan syntax).

## Usage

``` r
graph2lavaan(graph, nodes = V(graph)$name, ...)
```

## Arguments

- graph:

  A graph as an igraph object.

- nodes:

  Subset of nodes to be included in the model. By default, all the input
  graph nodes will be included in the output model.

- ...:

  Currently ignored.

## Value

A model in lavaan syntax.

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
# Graph (igraph object) to structural model in lavaan syntax
model <- graph2lavaan(sachs$graph)
cat(model, "\n")
#> zAkt~zPIP3
#> zAkt~zPKA
#> zErk~zMek
#> zErk~zPKA
#> zJnk~zPKA
#> zJnk~zPKC
#> zMek~zPKA
#> zMek~zPKC
#> zMek~zRaf
#> zP38~zPKA
#> zP38~zPKC
#> zPIP2~zPIP3
#> zPIP2~zPlcg
#> zPKC~zPIP2
#> zPKC~zPlcg
#> zPlcg~zPIP3
#> zRaf~zPKA
#> zRaf~zPKC
#> zPKA~~zPIP3 
```
