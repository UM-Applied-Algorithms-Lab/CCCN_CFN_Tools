% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KeyFunctions.R
\name{process_ptms_data}
\alias{process_ptms_data}
\title{Process PTMs Data}
\usage{
process_ptms_data(eu.sp.sed.ptms, sed.ptms.peps, AlldataPTMs_cor)
}
\arguments{
\item{eu.sp.sed.ptms}{A list of all PTMs.}

\item{sed.ptms.peps}{A list of SED PTMs peptides.}

\item{AlldataPTMs_cor}{A correlation matrix for all PTMs.}
}
\value{
A data frame containing PTMs gene correlation edges.
}
\description{
This function processes PTMs data, creates correlation networks, and constructs adjacency matrices.
}
\examples{
process_ptms_data(eu.sp.sed.ptms, sed.ptms.peps, AlldataPTMs_cor)
}
