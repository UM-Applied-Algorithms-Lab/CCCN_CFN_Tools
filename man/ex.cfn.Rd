% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_documentation.R
\docType{data}
\name{ex.cfn}
\alias{ex.cfn}
\title{Cluster Filtered Network Example}
\format{
Dataframe of 71 rows and 4 columns:
\describe{
\item{Gene.1}{First gene interacting}
\item{Gene.2}{Second gene interacting}
\item{Interaction}{The type of interaction}
\item{PPI.weight}{The weight of the interaction between these two genes as described by STRINGdb and GeneMANIA}
}
}
\source{
"created by MakeDBInput"
}
\usage{
ex.cfn
}
\description{
Subset of the PPI network in the cocluster correlation network.
}
\examples{
head(ex.cfn)
}
\keyword{internal}
