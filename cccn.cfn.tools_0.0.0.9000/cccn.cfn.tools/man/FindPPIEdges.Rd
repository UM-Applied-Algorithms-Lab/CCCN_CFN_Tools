% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindPPIEdges.R
\name{FindPPIEdges}
\alias{FindPPIEdges}
\title{Find PPI Edges}
\usage{
FindPPIEdges(
  cccn_matrix,
  db_filepaths = c(),
  gm.network = NA,
  ppi.network.name = "ppi.network"
)
}
\arguments{
\item{cccn_matrix}{dataframe of dataframes that represent the common clusters from the three distance calculations' clusters}

\item{db_filepaths}{vector of filepaths to data from other databases; defaults to empty vector}
}
\value{
A data frame of combined edges from STRINGdb and provided database entries
}
\description{
This function finds protein-protein interaction edges by combining STRINGdb and GeneMANIA databases.
}
\examples{
gmfile <- system.file("genemania", "genemania-interactions.txt", package = "cccn.cfn.tools", mustWork = TRUE)
cccn.cfn.tools:::ex.FindPPIEdges(ex.cccn_matrix, db_filepaths = c(), gm.network = NA, ppi.network.name = ex.ppi_network)
}
