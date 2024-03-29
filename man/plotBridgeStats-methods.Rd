% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllPlots.R
\docType{methods}
\name{plotBridgeStats,Bridge-method}
\alias{plotBridgeStats,Bridge-method}
\alias{plotBridgeStats}
\title{plotBridgeStats}
\usage{
\S4method{plotBridgeStats}{Bridge}(gbr, thr = 0.01)
}
\arguments{
\item{gbr}{A preprocessed \linkS4class{Bridge} class object.}

\item{thr}{A pvalue threshold.}
}
\value{
A ggplot object.
}
\description{
The \code{plotBridgeStats} generated a ggplot2 object for 
statistics inferred by the Bridge algorithm.
}
\examples{
# Load datasets used for demonstration
data(ogdata)

# Create an object of class 'Bridge' for H. sapiens
gbr <- newBridge(ogdata, phyloTree, refsp="9606", ogids=ogids)

# Run the Bridge algorithm
gbr <- runBridge(gbr)
gbr <- runPermutation(gbr, nPermutations=100)
plotBridgeStats(gbr)

}
