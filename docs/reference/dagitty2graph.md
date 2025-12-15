# Graph conversion from dagitty to igraph

Convert a dagitty object to a igraph object.

## Usage

``` r
dagitty2graph(dagi, verbose = FALSE, ...)
```

## Arguments

- dagi:

  A graph as a dagitty object ("dag" or "pdag").

- verbose:

  A logical value. If TRUE, the output graph is shown. This argument is
  FALSE by default.

- ...:

  Currently ignored.

## Value

An igraph object.

## Author

Mario Grassi <mario.grassi@unipv.it>

## Examples

``` r
# Conversion from igraph to dagitty  (and viceversa)
dagi <- graph2dagitty(sachs$graph, verbose = TRUE)

graph <- dagitty2graph(dagi, verbose = TRUE)

```
