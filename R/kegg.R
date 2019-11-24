# SEMgraph 0.3.3
#
# This is the SEMgraph package for statistical graph analysis using
# Structural Equation Models
# Typically, SEMgraph requires three input elements:
# 1) Interactome (e.g. KEGG singaling pathways or STRING PPI)
# 2) Quantitative data (e.g. GWAS, DNA-methylation arrays, RNA-seq)
# 3) Case/Control identifiers
#

#' @title KEGG signaling pathways collection
#'
#' @description Collection of KEGG signaling pathways stored as igraph objects.
#' @name kegg.pathways
#' @usage kegg.pathways
#' @docType data
#' @format
#' "kegg.pathways" is a collection of 214 KEGG signaling pathways, saved as graphNEL objects.
#' @source \url{https://www.kegg.jp/kegg}
#' @references
#' Kanehisa M1, Goto S (1999). KEGG: kyoto encyclopedia of genes and genomes. Nucleic Acid Research 28(1): 27-30. https://doi.org/10.1093/nar/27.1.29
#'
#' @examples
#' # Number of nodes per pathway
#' kegg.nodes <- unlist(lapply(kegg.pathways, graph::numNodes))
#' # Number of edges per pathway
#' kegg.edges <- unlist(lapply(kegg.pathways, graph::numEdges))
#' # Gene list per pathway
#' kegg.genes <- unlist(lapply(kegg.pathways, graph::nodes))
#' quantile(kegg.nodes)
#' quantile(kegg.edges)
#' length(unique(kegg.genes))  # Number of unique genes within the dataset
#' # Loading breast cancer KEGG network
#' i <- which(names(kegg.pathways) == "hsa05224_Breast_cancer")
#' bc.graph <- kegg.pathways[[i]]
#' bc.graph
NULL
