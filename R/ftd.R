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


#' @title Gene expression profiling of Frontotemporal Dementia (FTLD-U)
#' patients
#'
#' @description Microarray gene expression profiling of postmortem
#' brain samples from normal controls and FTLD-U patients with
#' progranulin gene mutations. Samples were taken from frontal cortex,
#' hippocampus, and cerebellum.
#' @name FTLDu_GSE13162
#' @usage FTLDu_GSE13162
#' @docType data
#' @format
#' "FTLDu_GSE13162" is a matrix of 12432 gene entrez IDs (rows) and 32
#' subjects (columns). Subjects are subdivided as follows:
#' \itemize{
#' \item rows 1 to 17, healthy subjects;
#' \item rows 18 to 32, FTLD-U subjects with progranulin mutations.
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13162}
#' @references
#' Chen-Plotkin AS1, Geser F, Plotkin JB, Clark CM, Kwong LK, Yuan W,
#' Grossman M, Van Deerlin VM, Trojanowski JQ, Lee VM. Variations in
#' the progranulin gene affect global gene expression in frontotemporal
#' lobar degeneration. Hum Mol Genet 2008 May 15;17(10):1349-62.
#' https://doi.org/10.1093/hmg/ddn023
#'
#' @examples
#' group <- c(rep(0, 17), rep(1, 15))
#' data <- FTLDu_GSE13162
#' # Set parameters and draw box plots
#' idx <- as.factor(group)
#' labels <- c("Control", "FTD")
#' pal <- c("red", "blue")
#' title <- "Gene Expression Intensities"
#' par(mar = c(3 + round(max(nchar(colnames(data)))/2), 4, 3, 1))
#' boxplot(data, boxwex = 0.6, notch = FALSE, outline = FALSE, las = 2,
#'         col = pal[idx], main = title)
#' legend("topleft", legend = labels, fill = pal, bty = "n", cex = 1)
NULL
