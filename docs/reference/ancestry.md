# Node ancestry utilities

Get ancestry for a collection of nodes in a graph. These functions are
wrappers for the original `SEMID` R package.

## Usage

``` r
ancestors(g, nodes)

descendants(g, nodes)

parents(g, nodes)

siblings(g, nodes)
```

## Arguments

- g:

  An igraph object.

- nodes:

  the nodes in the graph of which to get the ancestry.

## Value

a sorted vector of nodes.

## References

Rina Foygel Barber, Mathias Drton and Luca Weihs (2019). SEMID:
Identifiability of Linear Structural Equation Models. R package version
0.3.2. \<https://CRAN.R-project.org/package=SEMID/\>

## Examples

``` r
# Get all ancestors
an <- V(sachs$graph)[ancestors(sachs$graph, "Erk")]; an
#> + 8/11 vertices, named, from 39fcb44:
#> [1] PKA  PKC  PIP3 Mek  Raf  PIP2 Plcg Erk 

# Get parents
pa <- V(sachs$graph)[parents(sachs$graph, "PKC")]; pa
#> + 3/11 vertices, named, from 39fcb44:
#> [1] PKC  PIP2 Plcg

# Get descendants
de <- V(sachs$graph)[descendants(sachs$graph, "PKA")]; de
#> + 11/11 vertices, named, from 39fcb44:
#>  [1] PKA  PKC  PIP3 Mek  Raf  PIP2 Plcg Jnk  P38  Akt  Erk 

# Get siblings
sib <- V(sachs$graph)[siblings(sachs$graph, "PIP3")]; sib
#> + 5/11 vertices, named, from 39fcb44:
#> [1] PKA  PIP3 PIP2 Plcg Akt 
```
