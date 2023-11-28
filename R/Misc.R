
################################################################################
### Package documentation
################################################################################
#' @keywords internal
#' @title GeneBridge: large-scale evolutionary analysis of orthologous groups
#'
#' @description
#' GeneBridge is designed to assess the distribution of orthologous groups in 
#' a given species tree. It implements the Bridge algorithm to determine the 
#' evolutionary root of genes in large-scale evolutionary analysis.
#'
#' @details
#'
#' \tabular{ll}{
#' Package: \tab GeneBridge\cr
#' Type: \tab Software\cr
#' License: \tab GPL-3\cr
#' Maintainer: \tab Mauro Castro \email{mauro.a.castro@@gmail.com}\cr
#' }
#'
#' @section Index:
#' \tabular{ll}{
#' \link{newBridge}:
#' \tab Constructor of Bridge-class objects.\cr
#' \link{runBridge}:
#' \tab Run gene root inference with the Bridge algorithm.\cr
#' \link{getBridge}:
#' \tab Accessors for fetching slots from a Bridge object.\cr
#' }
#' Further information is available in the vignettes by typing
#' \code{vignette('GeneBridge')}. Documented topics are also available in
#' HTML by typing \code{help.start()} and selecting the GeneBridge package
#' from the menu.
#'
#'
"_PACKAGE"
#> [1] '_PACKAGE'

################################################################################
### Documentation for 'cogData'
################################################################################
#' @title A pre-processed dataset for the GeneBridge package
#'
#' @description Dataset used to demonstrate GeneBridge main functions.
#'
#' @format A data frame containing orthology annotation.
#'
#' @usage data(cogData)
#'
#' @source This package.
#' 
#' @details
#' 
#' The dataset consists of 4 R objects to be used for demonstration purposes 
#' in the geneplast vignette.
#' 
#' \describe{
#' 
#' \item{cogdata}{
#' A data frame with three columns listing orthology 
#' annotation retrieved from the STRING database (http://string-db.org/), 
#' release 9.1. Column 1 = Ensembl protein ID; column 2 = NCBI species ID; 
#' column 3 = OG ID. Note: This dataset is to be used for demonstration 
#' purposes only as it represents a subset of the STRING database; in order 
#' to reduce the dataset size, orthology annotation was mapped to genome 
#' stability genes (Castro et al.).
#' }
#' \item{sspids}{
#' A data frame containing the species listed in STRING database 
#' (http://string-db.org/), release 9.1. Column 1 = NCBI species ID; 
#' column 2 = NCBI species name;column 3 = species domain (eukaryotes).
#' }
#' \item{cogids}{
#' A one-column data.frame listing orthologous groups (OGs) available in 
#' the 'cogdata' object.
#' }
#' \item{phyloTree}{
#' An phylo-class object for the eukaryotes listed in the STRING database, 
#' release 9.1.
#' }
#' }
#' 
#' @references
#' Franceschini et al. STRING v9.1: protein-protein interaction networks, 
#' with increased coverage and integration. Nucleic Acids Research, 
#' 41(Database issue):D808-15, 2013. doi:10.1093/nar/gks1094.
#' 
#' @docType data
#' @keywords cogData
#' @name cogData
#' @aliases cogData
#' @return A pre-processed dataset.
#' @examples
#' data(cogData)
NULL

