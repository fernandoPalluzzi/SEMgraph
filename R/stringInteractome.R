# SEMgraph 0.3.3
#
# This is the SEMgraph package for statistical graph analysis using
# Structural Equation Models
# Typically, SEMgraph requires three input elements:
# 1) Interactome (e.g. KEGG singaling pathways or STRING PPI)
# 2) Quantitative data (e.g. GWAS, DNA-methylation arrays, RNA-seq)
# 3) Case/Control identifiers
#

#' @title STRING interactome
#'
#' @description STRING interactome version 10.5.
#' @name string
#' @usage string
#' @docType data
#' @format
#' "string" is an igraph network object of 9725 nodes and 170987 edges corresponding to the STRING interactome (version 10.5).
#' @source \url{https://string-db.org}
#' @references
#' Szklarczyk D, Gable AL, Lyon D, Junge A, Wyder S, Huerta-Cepas J, Simonovic M, Doncheva NT, Morris JH, Bork P, Jensen LJ, von Mering C (2019). STRING v11: protein-protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets. Nucleic Acids Res., 47: D607-613. https://doi.org/10.1093/nar/gky1131
#'
#' @examples
#' string
#' # STRING degrees of freedom
#' vcount(string)*(vcount(string) - 1)/2 - ecount(string)
NULL
