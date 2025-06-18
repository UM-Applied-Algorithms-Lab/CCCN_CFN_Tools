#' Bioplanet Example
#'
#' (description)
#'
#' @format Dataframe with 100 rows and 4 columns:
#' \describe{
#'   \item{}{}
#'   ...
#' }
#' @source source info (or `"manually created"`)
#' @examples
#' data(ex.bioplanet)
#' head(ex.bioplanet)
"ex.bioplanet"

#' CoCluster Correlation Network
#'
#' Network of relationships between genes based on clustering of ptms from the three distance metrics.
#'
#' @format Square matrix with 19 rows and columns:
#' \describe{
#'   \item{rows/columns}{Gene}
#' }
#' @source "created by MakeCorrelationNetwork"
#' @examples
#' data(ex.cccn_matrix)
#' head(ex.cccn_matrix)
"ex.cccn_matrix"

#' Pathways List Example
#'
#' (description)
#'
#' @format List of 1657 elements:
#' \describe{
#'   \item{pathway}{A pathway name with the list of genes involved}
#' }
#' @source source info (or `"manually created"`)
#' @examples
#' data(ex.pathways.list)
#' head(ex.pathways.list)
"ex.pathways.list"

#' Full GeneMANIA Edges Example
#'
#' Network of interactions from the GeneMANIA database.
#'
#' @format Dataframe of 601,552 rows and 10 columns:
#' \describe{
#'   \item{data.type}{Type of interaction; typically co-expression}
#'   \item{highlight}{UNKNOWN}
#'   \item{interaction}{UNKNOWN}
#'   \item{name}{The names of the two genes interacting as well as the data type}
#'   \item{networks}{Source of information}
#'   \item{normalized.max.weight}{UNKNOWN}
#'   \item{raw.weights}{A number weighting the strength of the interaction}
#'   \item{selected}{UNKNOWN}
#'   \item{shared.interaction}{UNKNOWN}
#'   \item{shared.name}{UNKNOWN}
#'   ...
#' }
#' @source \url{https://apps.cytoscape.org/apps/genemania}
#' @examples
#' data(ex.gm.edges)
#' head(ex.gm.edges)
"ex.gm.edges"

#' Full GeneMANIA Network Example
#'
#' Network of relationships between genes based on information from the GeneMANIA database.
#'
#' @format Dataframe of 601,552 rows and 3 columns:
#' \describe{
#'   \item{}{}
#'   ...
#' }
#' @source \url{https://apps.cytoscape.org/apps/genemania}
#' @examples
#' data(ex.gm.network)
#' head(ex.gm.network)
"ex.gm.network"

#' GeneMANIA Nodes Example
#'
#' Dataframe from GeneMANIA relating GeneMANIA IDs to common gene names.
#'
#' @format Dataframe of 2762 rows and 16 columns:
#' \describe{
#'   \item{}{}
#'   ...
#' }
#' @source \url{https://apps.cytoscape.org/apps/genemania}
#' @examples
#' data(ex.gm.nodes)
#' head(ex.gm.nodes)
"ex.gm.nodes"

#' GeneMANIA Empty Edges Example
#'
#' Empty GeneMANIA structure. TODO: fix.
#'
#' @format Dataframe of 0 rows and 4 columns:
#' \describe{
#'   \item{}{}
#'   ...
#' }
#' @source \url{https://genemania.org}
#' @examples
#' data(ex.gm_edges)
#' head(ex.gm_edges)
"ex.gm_edges"

#' List of Common Clusters Example
#'
#' (description)
#'
#' @format List of 12 variables:
#' \describe{
#'   \item{}{}
#'   ...
#' }
#' @source source info (or `"manually created"`)
#' @examples
#' data(ex.list.common)
#' head(ex.list.common)
"ex.list.common"

#' Nodenames Example
#'
#' List of gene names in dataset.
#'
#' @format Dataframe of 10 rows and 1 column:
#' \describe{
#'   \item{}{}
#'   ...
#' }
#' @source source info (or `"manually created"`)
#' @examples
#' data(ex.nodenames)
#' head(ex.nodenames)
"ex.nodenames"

#' PPI Network Example
#'
#' Network of relationships between genes based on PPI database STRINGdb.
#'
#' @format Dataframe of 1 row and 4 columns:
#' \describe{
#'   \item{}{}
#'   ...
#' }
#' @source source info (or `"manually created"`)
#' @examples
#' data(ex.ppi_network)
#' head(ex.ppi_network)
"ex.ppi_network"

#' Correlation Matrix Example
#'
#' Matrix relating the correlation values between ptms based on Spearman Similarity.
#'
#' @format Square matrix of 100 rows and columns:
#' \describe{
#'   \item{}{}
#'   ...
#' }
#' @source source info (or `"manually created"`)
#' @examples
#' data(ex.ptm.correlation.matrix)
#' head(ex.ptm.correlation.matrix)
"ex.ptm.correlation.matrix"

#' PTM Table Example
#'
#' Post-translational modification data under certain environmental conditions.
#'
#' @format Dataframe of 100 rows and 70 columns:
#' \describe{
#'   \item{}{}
#'   ...
#' }
#' @source source info (or `"manually created"`)
#' @examples
#' data(ex.ptmtable)
#' head(ex.ptmtable)
"ex.ptmtable"

#' TSNE matrices Example
#'
#' List of of lists of dataframes of ptm positions in 3D space to show their relationships with each other based on distance.
#'
#' @format List of 3 lists of 6 (Euclidean and SED) or 88 (Spearman) dataframes of varying sizes:
#' \describe{
#'   \item{}{}
#'   ...
#' }
#' @source source info (or `"manually created"`)
#' @examples
#' data(ex.tsne.matrices)
#' head(ex.tsne.matrices)
"ex.tsne.matrices"

