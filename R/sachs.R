#  SEMgraph library
#  Copyright (C) 2019-2021 Mario Grass; iFernando Palluzzi 
#  e-mail: <fernando.palluzzi@gmail.com>
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


#' @title Sachs multiparameter flow cytometry data and consensus model
#'
#' @description Flow cytometry data and causal model from Sachs et al. (2005).
#' @name sachs
#' @usage sachs
#' @docType data
#' @format
#' "sachs" is a list of 5 objects:
#' \enumerate{
#' \item "rawdata", a list of 14 data.frames containing raw flow cytometry
#' data (Sachs et al., 2005);
#' \item "graph", consensus signaling network;
#' \item "model", consensus model (lavaan syntax);
#' \item "pkc", data.frame of 1766 samples and 11 variables, containing
#' cd3cd28 (baseline) and pma (PKC activation) data;
#' \item "group", a binary group vector, where 0 is for cd3cd28 samples
#' (n = 853) and 1 is for pma samples (n = 913).
#' \item "details", a data.frame containing dataset information.
#' }
#' @source \doi{10.1126/science.1105809}
#'
#' @references
#'
#' Sachs K, Perez O, Pe'er D, Lauffenburger DA, Nolan GP (2019).
#' Causal Protein-Signaling Networks Derived from Multiparameter
#' Single-Cell Data. Science, 308(5721): 523-529.
#'
#' @examples
#' # Dataset content
#' names(sachs$rawdata)
#' dim(sachs$pkc)
#' table(sachs$group)
#' cat(sachs$model)
#' gplot(sachs$graph)
#'
NULL
