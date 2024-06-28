#---------------------------------------------------------#
#  Loading KEGG and Reactome pathway data using graphite  #
#---------------------------------------------------------#


# Loading graphite and preparing the environment

library(graphite)
library(org.Hs.eg.db)          # Homo sapiens annotations
library(org.Mm.eg.db)          # Mus Musculus annotations
library(igraph)
library(SEMgraph)


#' @title Import pathways and generate a reference network
#'
#' @description Utility to create pathway lists as igraph objects and 
#' interaction networks from Reactome, KEGG, and other pathway databases.
#'
#' @param graph An igraph object.
#' @param ... Currently ignored.
#'
#' @details This function uses \code{graphite} to download and preprocess 
#' network data from pathway databases. The output is then created using 
#' igraph and SEMgraph utilities.
#'
#' @return A list of 2 objects:
#' \enumerate{
#' \item a list of pathways ad igraph objects;
#' \item the union of graphs in the pathway list.
#' }
#'
#' @seealso \code{\link[graphite]{pathways}} and 
#' \code{\link[graphite]{pathwayDatabases}}
#'
#' @import igraph
#' @import graphite
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#'
#' @export
#'
#' @references
#'
#' Sales G, Calura E, Cavalieri D, Romualdi C (2012). graphite - a Bioconductor 
#' package to convert pathway topology to gene network. 
#' BMC Bioinformatics. 
#' <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-20>
#'
#' Csardi G, Nepusz T (2006). The igraph software package for complex network 
#' research. InterJournal, Complex Systems, 1695.
#' <https://igraph.org>
#'
#' Grassi M, Palluzzi F, Tarantino B (2022). SEMgraph: An R Package for Causal Network
#' Analysis of High-Throughput Data with Structural Equation Models.
#' Bioinformatics, 38 (20), 4829–4830 <https://doi.org/10.1093/bioinformatics/btac567>
#'
#' @examples
#' 
#' # Create KEGG reference pathway list and reference network for Homo sapiens
#' kegg.hs <- loadPathwayData("hsapiens", "kegg", id_type = "ENTREZID")
#' 
#' # Inspect results
#' length(kegg.hs$pathways)
#' kegg.hs$network
#'
loadPathwayData <- function(organism, db, id_type = "ENTREZID", lcc = TRUE) {
  
  # Import pathway from database
  pathways <- pathways(organism, db)
  pathway.list <- list()
  n <- length(pathways)
  
  for (i in 1:n) {
    
    pw.name <- names(pathways)[[i]]
    message(paste0("Pathway ", i, " of ", n, ": ", pw.name))
    
    # Node ID conversion
    pw <- graphite::convertIdentifiers(pathways[[i]], id_type)
    
    # Conversion pathway -> GraphNEL -> igraph
    G <- igraph::igraph.from.graphNEL(pathwayGraph(pw))
    
    # Removing node name prefix
    V(G)$name <- sub(paste0(id_type, ":"), "", V(G)$name)
    
    # Adding graph to the list
    pathway.list[[i]] <- G
    names(pathway.list)[[i]] <- pw.name
  }
  
  # Generating reference network
  reference <- igraph::graph.union(pathway.list)
  if (lcc) {
    reference <- SEMgraph::properties(reference)[[1]]
  }
  return(list(pathways = pathway.list, network = reference))
}

