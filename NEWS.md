## Version 1.2.4 Release Notes
* Restored `kegg.RData` and `kegg.pathways.RData` (November 2021) to ensure
retro‑compatibility with benchmarked examples.

* Added a new `loadPathways()` function to retrieve the latest pathway lists  
and their union (interactome) as *igraph* objects, using a wrapper for the
*graphite* package from the Bioconductor project.

* Added a tutorial vignette "Get Started", that was only on the Github before.

* Various fixed bugs discovered after the release 1.2.3.

## Version 1.2.3 Release Notes
* Update `kegg.RData` and `kegg.pathways.RData` (February 2025).

* Various fixed bugs discovered after the release 1.2.2.

## Version 1.2.2 Release Notes
* Delete `predictSink()` function. A general function for SEM-based
out-of-sample prediction is now included in the *SEMdeep* package,
which uses Deep Neural Network (DNN) and Machine Learning(ML) algorithms,
has been released on CRAN: 10.32614/CRAN.package.SEMdeep

* Various fixed bugs discovered after the release 1.2.1.

## Version 1.2.1 Release Notes
* Added new `predictSink()` function for SEM-based out-of-sample prediction
of (observed) response y-variables (sink nodes) given the values
of (observed) x-variables (source and mediator) nodes from the fitted
graph structure.

* Added new `transformData()` function implementing various data trasformation
methods to perform optimal scaling for ordinal or nominal data, and to help
relax the assumption of normality (gaussianity) for continuous data. 

* Update `kegg.RData` and `kegg.pathways.RData` (February 2024).

* Various fixed bugs discovered after the release 1.2.0.

## Version 1.2.0 Release Notes
* Version 1.2.0 is a major release with several new features, including:

* `SEMrun()` function.  The algo ="cggm" based on high-dimensional GGGM is now
implemented with the de-sparsified (de-biased) nodewise LASSO procedure 
applied on a Gaussian DAG model. The overall indices "deviance/df" and "srmr"
are now computed using the observed correlation matrix also in p > n regime,
where the estimated parameters are computed using the "regularized" (lambda
corrected) correlation matrix.

* `SEMbap()` function. New deconfounding methods to adjust the data matrix
by removing latent sources of confounding encoded in them are implemented.
The selected methods are either based on: (i) Bow-free Acyclic Paths (BAP)
search (dalgo = "cggm" or "glpc"), (ii) LVs proxies as additional source
nodes of the data matrix, Y (dalgo = "pc" or "glpc") or (iii) spectral
transformation of Y (dalgo = "pc" or "trim").

* `SEMdag()` function. New two-step DAG estimation from an input (or empty) graph,
using in step 1) graph topological order or bottom-up search order, and in
step 2) parent recovery with the LASSO-based algorithm are implemented.
The estimate linear order are obtained from a priori graph topological vertex
(LO = "TO") or level (LO = "TL") ordering, or with a data-driven vertex or 
level Bottom-up (LO = "BU") based on "glasso" residual variance ordering.
The Top-Down (LO = "TD") is removed, being the BU more efficient to implement
the topological search order.

* `Shipley.test()` function. Added new argument cmax = Inf (default). This
parameter can be used to perform only those tests where the number of
conditioning variables does not exceed the given value. Output of the
data.frame "dsep" has the same format of the localCI.test() function.

* Various fixed bugs discovered after the release 1.1.3.

## Version 1.1.3 Release Notes
* Added in `SEMrun()` function the argumet SE = c("standard" or "none"), if
algo = "lavaan".

* Added in `SEMrun()` function the bootstrap resampling of SE (95% CI), and
new argoment n_rep = 1000 (default) to set the bootstrap samples or permutation
flip, if algo = "ricf".

* Added in `SEMrun()` function the de-sparsified SE (95% CI) of omega parameters 
(the elements of the precision matrix), if algo = "cggm".

* Added new `parameterEstimates()` function for parameter estimates output
of a fitted SEM for RICF and CGGM algorithms similar to lavaan.

* Updating `summary.RICF()` and `summary.GGM()` functions with `parameterEstimates()`.

* Various fixed bugs discovered after the release 1.1.2.

## Version 1.1.2 Release Notes
* Added new `SEMtree()` function for tree-based structure learning methods.
Four methods with graph (type= "ST" or "MST") and data-driven (type = "CAT"
or "CPDAG") algorithms are implemented.

* Deprecated `activeModule()` and `corr2graph()` functions in favor of new `SEMtree()`
function. 

* Added new `dagitty2graph()` function for conversion from a *dagitty* graph object
to an igraph object.

* Added new `localCI.test()` function for local conditional indipendence (CI)
test of missing edges from an acyclic graph. This function is a wrapper to
the `localTests()` function from package *dagitty*. 

* Added new arguments for `SEMace()` function: type = c("parents", "minimal",
"optimal") to choose the conditioning set Z of Y over X; effect = c("all",
"source2sink", "direct",) to choose the type of X to Y effect. 

* Added new argument for `SEMdci()` function: type = "ace" from `SEMace()` function
with fixed type="parents", and effect="direct".

* Change `mergeGraph()` function. Now the function combines groups of graph
nodes using hierarchical clustering with prototypes derived from protoclust
package or custom membership attribute (e.g., cluster membership derived from
`clusterGraph()` function).  

* Delete argument seed = c(0.05, 0.5, 0.5) in the function `weigthGraph()`. Now
if group is NOT NULL also node weighting is actived, and node weights correspond
to the sign and P-value of the z-test = b/SE(b) from glm(node ~ group).

* Various fixed bugs discovered after the release 1.1.0.

## Version 1.1.0 Release Notes
* Version 1.1.0 is a major release with significant changes:

* Added new arguments for `SEMdag()` function: LO = "TO" or "TD" for knowledge-based
topological order (TO) or data-driven top-down order (TD), and penalty = TRUE or
FALSE, binary penalty factors can be applied to each L1-coefficient.

* Deprecated `extendGraph()` in favor of new `resizeGraph()` function, that 
re-sized graph, removing edges or adding edges/nodes if they are present
or absent in a given reference network.

* Change `modelSerch()`, interactive procedure is out, and now a three step
procedure is implemented for search strategies with new `SEMdag()` and `resizeGraph()`
functions.

* Change `SEMgsa()` deleting D,A,E p-values with more performing activation and
inhibition pvalues.

* Added argument MCX2= TRUE or FALSE for `Shipley.test()` function, a Monte Carlo
P-value of the combined C test.

* Added new `SEMdci()` function for differentially connected genes inference.

* Change `properties()`, now extracted components are order by component sizes.

* Change argument q = q-quantile with q = 1-top/vcount(graph) in `activeModule()`
function, now the induced graph for the "rwr" and "hdi" algorithms is defined
by the top-n ranking nodes.

* Various fixed bugs.

## Version 1.0.3 Release Notes
* First stable version on CRAN.

* Update `kegg.RData` (November, 2021).

* Added `kegg.pathways.RData` (November, 2021).

* Added *pkgdown* website.

* Various fixed bugs.

## Version 1.0.0 Release Notes
* First stable version on GitHub
