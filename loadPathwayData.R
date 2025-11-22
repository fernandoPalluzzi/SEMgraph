#---------------------------------------#
#  Loading pathway data using graphite  #
#---------------------------------------#

#' @title Import pathways and generate a reference network
#'
#' @description Utility to create pathway lists as igraph objects and 
#' interaction networks from Reactome, KEGG, and other pathway databases.
#'
#' @param organism A string indicating the source organism. Please, check 
#' the \code{pathwayDatabases} function from graphite to list the available 
#' datasets. 
#' @param db String indicating the database name. Please, check the 
#' \code{pathwayDatabases} function from graphite to list the available 
#' datasets.
#' @param id_type Gene ID type. The default is set to "ENTREZID" (standard 
#' SEMlearn nomenclature). A common choice could be "SYMBOL", for 
#' official gene symbols.
#' @param lcc A logical value. If TRUE (default), the reference network 
#' will only include the largest connected component. It will include all 
#' disconnected components otherwise.
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
loadPathwayData0 <- function(db, organism = "hsapiens", id_type = "ENTREZID", lcc = TRUE, ...) {
  
  # Import pathway from database
  pathways <- pathways(organism, db)
  pathway.list <- list()
  n <- length(pathways)
  del <- NULL
  
  for (i in 1:n) { #i=51
    pw.name <- names(pathways)[[i]]
	#if (i %% 10 == 0) message("  Processed ", i, "/", n, " pathways")
    #message(paste0("Pathway ", i, " of ", n, ": ", pw.name))
    
    # Node ID conversion
    pw <- graphite::convertIdentifiers(pathways[[i]], id_type)
    
    # Conversion pathway -> GraphNEL -> igraph
    G <- igraph::graph_from_graphnel(pathwayGraph(pw))
	G <- simplify(G, remove.loops = TRUE)
    G <- G - E(G)[which_mutual(G)]
  	
	# Skip small graph
    if (igraph::vcount(G) <= 5 || igraph::ecount(G) == 0) {
	 message(paste0("Delete pathway ", i, " of ", n, ": ", pw.name))
	 del <- c(del, i); next
	}
		    
    #Removing node name prefix and edge weight
	V(G)$name <- sub(paste0(id_type, ":"), "", V(G)$name)
	G <- delete_edge_attr(G, "weight")
	
	# Adding graph to the list
    pathway.list[[i]] <- G
    names(pathway.list)[[i]] <- pw.name
  }
  
  # Generating reference network
  reference <-  Reduce(igraph::union, pathway.list)
  if (lcc) reference <- properties(reference)[[1]] #SEMgraph
  
  # Add node labels to all graphs
  add_labels <- function(g) {
    V(g)$label <- suppressMessages(mapIds(org.Hs.eg.db, V(g)$name, "SYMBOL", id_type))
    g
  }
  pathway.list <- lapply(pathway.list[-del], add_labels)
  reference <- add_labels(reference)
  
  return(list(pathways = pathway.list, network = reference))
}
