% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClusterFilteredNetwork.R
\name{ClusterFilteredNetwork}
\alias{ClusterFilteredNetwork}
\title{Filter PPI Network by Known Clusters}
\usage{
ClusterFilteredNetwork(cccn.matrix, ppi.network, cfn.name = "cfn")
}
\arguments{
\item{cccn.matrix}{Matrix representing the common clusters from the three distance calculations' clusters}

\item{ppi.network}{A data frame of combined edges from STRINGdb and provided database entries}

\item{cfn.name}{Desired name of the output cluster filtered network}
}
\value{
cluster filtered network
}
\description{
All edges in the PPI network are checked against the cocluster correlation network to ensure existence and a non-zero weight.
If either of these conditions are not met, then it will be removed from the list of PPI edges.
This new, cluster filtered network is then assigned to the global namespace.
}
\examples{
ClusterFilteredNetwork(ex.cccn.matrix, ex.ppi.network, cfn.name = "example.cfn")
}
