## Version 1.1.3 Release Notes

* Added new parameterEstimates() function for parameter estimates output
of a fitted SEM for RICF and CGGM algorithms similar to lavaan.

* Updating summary.RICF() and summary.GGM() functions with parameterEstimates()
(and hidden functions SEMricf(), SEMricf2(), SEMggm(), and SEMggm2(), ...)

* Added in SEMrun() function the bootstrap resampling of SE (95% CI), and
new argoment n_rep = 1000 (default) to set the bootstrap samples or permutation
flip, if algo = "ricf" and the argumet SE = c("standard" or "none"), if
algo = "lavaan".

* Fixed bugs in cplot, colorGraph, SEMtree, and others.


## Version 1.1.2 Release Notes

* Added new SEMtree() function for tree-based structure learning methods.
Four methods with graph (type= "ST" or "MST") and data-driven (type = "CAT"
or "CPDAG") algorithms are implemented.

* Deprecated activeModule() and corr2graph() functions in favor of new SEMtree()
function. 

* Added new dagitty2graph() function for conversion from a dagitty graph object
to an igraph object.

* Added new localCI.test() function for local conditional indipendence (CI)
test of missing edges from an acyclic graph. This function is a wrapper to
the function localTests() from package dagitty. 

* Added new arguments for SEMace() function: type = c("parents", "minimal",
"optimal") to choose the conditioning set Z of Y over X; effect = c("all",
"source2sink", "direct",) to choose the type of X to Y effect. 

* Added new argument for SEMdci() function: type = "ace" from ACE function()
with fixed type= "parents", and effect="direct"

* Change mergeGraph() function. Now the function combines groups of graph
nodes using hierarchical clustering with prototypes derived from protoclust
package or custom membership attribute (e.g., cluster membership derived from
function clusterGraph()).  

* Delete argument seed = c(0.05, 0.5, 0.5) in the function weigthGraph(). Now
if group is NOT NULL also node weighting is actived, and node weights correspond
to the sign and P-value of the z-test = b/SE(b) from glm(node ~ group).

* Various fixed bugs


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

* Added argument MCX2= TRUE or FALSE for Shipley.test() function, a Monte Carlo
P-value of the combined C test

* Added new SEMdci() function for differentially connected genes inference

* Change properties(), now extracted components are order by component sizes

* Change argument q = q-quantile with q = 1-top/vcount(graph) in activeModule()
function, now the induced graph for the "rwr" and "hdi" algorithms is defined
by the top-n ranking nodes.

* Various fixed bugs


## Version 1.0.3 Release Notes

* First stable version on CRAN

* Update kegg.RData (November, 2021)

* Added kegg.pathways.RData (November, 2021)

* Added pkgdown website

* Various fixed bugs


## Version 1.0.0 Release Notes

* First stable version on GitHub
