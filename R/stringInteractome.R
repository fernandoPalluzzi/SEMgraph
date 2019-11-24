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
