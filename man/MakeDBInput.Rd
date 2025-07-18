% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GatherPPIData.R
\name{MakeDBInput}
\alias{MakeDBInput}
\title{Make Database Input File}
\usage{
MakeDBInput(gene.cccn, file.path.name = "db_nodes.txt")
}
\arguments{
\item{gene.cccn}{A matrix showing strength of relationships between proteins using the common clusters between the three distance metrics (Euclidean, Spearman, and Combined (SED))}

\item{file.path.name}{Path for the output file; defaults to db_nodes.txt}
}
\value{
A file with all of the gene names which can be copy and pasted into the GeneMania cytoscape app, data frame of the names of the genes
}
\description{
This function outputs a file consisting entirely of gene names, each produced on a new line. This data can be copy and pasted into
a database input in order to get protein-protein interaction data.
}
\examples{
#MakeDBInput(ex.nodenames)
cat(ex.nodenames[[1]], sep = '\n')
}
