#  SEMgraph library
#  Copyright (C) 2019 Fernando Palluzzi; Mario Grassi
#  e-mail: <fernando.palluzzi@gmail.com>
#  University of Pavia, Department of Brain and Behavioral Sciences
#  Via Bassi 21, Pavia, 27100 Italy

#  SEMgraph is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  SEMgraph is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------- #


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
