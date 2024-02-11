#  SEMgraph library
#  Copyright (C) 2019-2024 Mario Grassi; Fernando Palluzzi; Barbara Tarantino
#  e-mail: <mario.grassi@unipv.it>
#  University of Pavia, Department of Brain and Behavioral Sciences
#  Via Bassi 21, 27100 Pavia, Italy

#  SEMgraph is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  SEMgraph is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

# -------------------------------------------------------------------- #

#' @title KEGG pathways
#'
#' @description KEGG pathways extracted using the \code{ROntoTools}
#' R package (update: February, 2024).
#' @name kegg.pathways
#' @usage kegg.pathways
#' @docType data
#' @format
#' "kegg.pathways" is a list of 227 igraph objects corresponding to the
#' KEGG pathways.
#' @source \url{https://www.genome.jp/kegg/}
#' @references
#' 
#' Kanehisa M, Goto S (1999). KEGG: kyoto encyclopedia of genes and genomes. 
#' Nucleic Acid Research 28(1): 27-30. 
#' <https://doi.org/10.1093/nar/27.1.29>
#'
#' Calin Voichita, Sahar Ansari and Sorin Draghici (2023).
#' ROntoTools: R Onto-Tools suite. R package version 2.30.0.
#' 
#' @examples
#' 
#' \donttest{
#' library(igraph)
#' 
#' # KEGG pathways
#' names(kegg.pathways)
#' 
#' i<-which(names(kegg.pathways)=="Type II diabetes mellitus");i
#' ig<- kegg.pathways[[i]]
#' summary(ig)
#' V(ig)$name
#' E(ig)$weight
#' 
#' gplot(ig, l="fdp", psize=50, main=names(kegg.pathways[i]))
#' 
#' }
#'
NULL
