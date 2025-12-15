# Graph conversion from igraph to dagitty

Convert an igraph object to a dagitty object.

## Usage

``` r
graph2dagitty(graph, graphType = "dag", verbose = FALSE, ...)
```

## Arguments

- graph:

  A graph as an igraph or as an adjacency matrix.

- graphType:

  character, is one of "dag" (default)' or "pdag". DAG can contain the
  directed (-\>) and bi-directed (\<-\>) edges, while PDAG can contain
  the edges: -\>, \<-\>, and the undirected edges (–) that represent
  edges whose direction is not known.

- verbose:

  A logical value. If TRUE, the output graph is shown. This argument is
  FALSE by default.

- ...:

  Currently ignored.

## Value

A dagitty object.

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
# Graph as an igraph object to dagitty object
G <- graph2dagitty(sachs$graph)
plot(dagitty::graphLayout(G))

```
