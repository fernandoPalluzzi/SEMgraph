## Version 1.1.0 Release Notes
* Major release with significant changes:

* Added new arguments for SEMdag() function: LO = "TO" or "TD" for knowledge-based
topological order or data-driven top-down order, and penalty = TRUE or FALSE,
separate penalty factors can be applied to each coefficient

* Deprecated extendGraph() in favor of new resizeGraph() function, that 
re-sized graph, removing edges or adding edges/nodes if they are present
or absent in a given reference network

* Change modelSerch(), interactive procedure is out, and now a three step
procedure is implemented for search strategies with new SEMdag() and resizeGraph
functions 

* Change SEMgsa() deleting D,A,E p-values with more performing activation and
inhibition pvalues

* Added argument MCX2 = TRUE or FALSE for Shipley.test() function, a Monte Carlo
P-value of the combined C test

* Added new SEMdci() function for differentially connected genes inference

* Change properties(), now extracted components are order by component sizes

* Change argument q = q-quantile with q = 1-top/vcount(graph) in activeModule()
function, now the induced graph for the "rwr" and "hdi" algorithms is defined
by the top-n ranking nodes.

* Various fixed bugs

## Version 1.0.5 Release Notes
* Fixed bugs for an update from 2.0.0 to 2.1.0 of GGMncv

## Version 1.0.4 Release Notes
* Fixed bugs for CRAN checking

## Version 1.0.3 Release Notes
* Fixed bugs for CRAN checking
* Update kegg.RData (November, 2021)
* Added kegg.pathways.RData (November, 2021)
* Added pkgdown website

## Version 1.0.2 Release Notes
* Change `orientEdges` function
* Added `graph2dagitty` function
* Fixed bugs

## Version 1.0.1 Release Notes
* Change parallel computation using `pbapply` package
* Fixed bugs 

## Version 1.0.0 Release Notes
* First stable version in GitHub
