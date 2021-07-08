#  SEMgraph library
#  Copyright (C) 2019-2021 Fernando Palluzzi; Mario Grassi
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

#' @title Amyotrophic Lateral Sclerosis (ALS) dataset
#'
#' @description Expression profiling through high-throughput sequencing 
#' (RNA-seq) of 139 ALS patients and 21 healthy controls (HCs), 
#' from Tam et al. (2019).
#' @name alsData
#' @usage alsData
#' @docType data
#' @format
#' alsData is a list of 4 objects:
#' \enumerate{
#' \item "graph", ALS graph as the largest connected component of the 
#' "Amyotrophic lateral sclerosis (ALS)" pathway from KEGG database;
#' \item "exprs", a matrix of 160 rows (subjects) and 303 columns (genes)
#' extracted from the original 17695. This subset includes genes from 
#' KEGG ALS signaling pathway, MAPK signaling pathway, and Steroid 
#' biosynthesis pathway, needed to run SEMgraph examples.
#' Raw data from the GEO dataset GSE124439 (Tam et al., 2019) were 
#' pre-processed applying batch effect correction, using the sva R package 
#' (Leek et al., 2012), to remove data production center and brain area 
#' biases. Using multidimensional scaling-based clustering, ALS-specific 
#' and an HC-specific clusters were generated. Misclassified samples were 
#' blacklisted and removed from the current dataset;
#' \item "group", a binary group vector of 139 ALS subjects (1) and 21 
#' healthy controls (0);
#' \item "details", a data.frame reporting information about included 
#' and blacklisted samples.
#' }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124439/}
#'
#' @references
#'
#' Tam OH, Rozhkov NV, Shaw R, Kim D et al. (2019). Postmortem Cortex 
#' Samples Identify Distinct Molecular Subtypes of ALS: Retrotransposon 
#' Activation, Oxidative Stress, and Activated Glia. Cell Repprts, 
#' 29(5):1164-1177.e5. <https://doi.org/10.1016/j.celrep.2019.09.066>
#' 
#' Jeffrey T. Leek, W. Evan Johnson, Hilary S. Parker, Andrew E. Jaffe, 
#' and John D. Storey (2012). The sva package for removing batch effects 
#' and other unwanted variation in high-throughput experiments. 
#' Bioinformatics. Mar 15; 28(6): 882-883. 
#' <https://doi.org/10.1093/bioinformatics/bts034>
#'
#' @examples
#' alsData$graph
#' dim(alsData$exprs)
#' table(alsData$group)
#'
NULL
