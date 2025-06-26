#' Bioplanet Example
#'
#' (description)
#'
#' @keywords internal
#'
#' @format Dataframe with 100 rows and 4 columns:
#' \describe{
#'   \item{PATHWAY_ID}{ID of the pathway in bioplanet}
#'   \item{PATHWAY_NAME}{Common name of the pathway}
#'   \item{GENE_ID}{ID of the gene in bioplanet}
#'   \item{GENE_SYMBOL}{Common name of the gene}
#' }
#' @source \url{https://tripod.nih.gov/bioplanet/}
#' @examples
#' data(ex.bioplanet)
#' head(ex.bioplanet)
"ex.bioplanet"



# DO NOT TOUCH THIS. DO NOT. IDK WHY. BUT DON'T TOUCH. DO NOT MAKE INTERNAL.
# TODO: fix

#' CoCluster Correlation Network
#'
#' Network of relationships between genes based on clustering of ptms from the three distance metrics.
#'
#' @format Square matrix with 19 rows and columns:
#' \describe{
#'   \item{rows/columns}{Gene}
#'   \item{values}{Strength of correlation}
#' }
#' @source "created by MakeCorrelationNetwork"
#' @examples
#' data(ex.cccn.matrix)
#' head(ex.cccn.matrix)
"ex.cccn.matrix"



#' Pathways List Example
#'
#' (description)
#'
#' @keywords internal
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
#' @keywords internal
#'
#' @format Dataframe of 601,552 rows and 10 columns:
#' \describe{
#'   \item{data.type}{Type of interaction; typically co-expression}
#'   \item{highlight}{Internal flag}
#'   \item{interaction}{Type of interaction}
#'   \item{name}{The names of the two genes interacting as well as the data type}
#'   \item{networks}{Source of information}
#'   \item{normalized.max.weight}{raw.weights scaled to be consistent}
#'   \item{raw.weights}{A number weighting the strength of the interaction}
#'   \item{selected}{Internal flag}
#'   \item{shared.interaction}{Type of interaction}
#'   \item{shared.name}{The names of the two genes interacting as well as the data type}
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
#' @keywords internal
#'
#' @format Dataframe of 601,552 rows and 3 columns:
#' \describe{
#'   \item{Gene.1}{The gene interacting with Gene.2}
#'   \item{Gene.2}{The gene interacting with Gene.1}
#'   \item{raw.weights}{Strength of the interaction}
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
#' @keywords internal
#'
#' @format Dataframe of 2762 rows and 16 columns:
#' \describe{
#'   \item{annotation.name}{Note about the gene}
#'   \item{annotations}{Processes associated with this gene according to Gene Ontology}
#'   \item{Ensemble.Gene.ID}{ID of the gene in Ensemble}
#'   \item{Entrez.Gene.ID}{ID of the gene in Entrez}
#'   \item{gene.name}{Common name fo the gene}
#'   \item{log.score}{Association strength to the set of entered genes on a log scale}
#'   \item{name}{Name of the gene in GeneMANIA}
#'   \item{node.type}{Why the node is in the list}
#'   \item{query.term}{What was searched for; typically the common name of the gene}
#'   \item{RefSeq.Protein.ID}{ID of the gene in RefSeq}
#'   \item{score}{Predicted association strength to the set of entered genes}
#'   \item{selected}{Internal flag}
#'   \item{shared.name}{Name of the gene in GeneMANIA}
#'   \item{Synonym}{Gene synonym in another database}
#'   \item{Uniprot.ID}{ID of the gene in Uniprot}
#' }
#' @source \url{https://apps.cytoscape.org/apps/genemania}
#' @examples
#' data(ex.gm.nodes)
#' head(ex.gm.nodes)
"ex.gm.nodes"


# DO NOT TOUCH I DON'T TRUST ITS ROXYGEN
# NO PROOF YET
#' GeneMANIA Empty Edges Example
#'
#' Empty GeneMANIA structure. TODO: fix.
#'
#' @format Dataframe of 0 rows and 4 columns:
#' \describe{
#'   \item{Gene.1}{The gene interacting with Gene.2}
#'   \item{Gene.2}{The gene interacting with Gene.1}
#'   \item{Weight}{Strength of the interaction}
#'   \item{Network.group}{Type of interaction}
#' }
#' @source \url{https://genemania.org}
#' @examples
#' data(ex.gm_edges)
#' head(ex.gm_edges)
"ex.gm_edges"


# DO NOT TOUCH THIS ALSO BROKE IT
# TODO: FIGURE OUT ??????
#' List of Common Clusters Example
#'
#' (description)
#'
#' @format List of 12 lists:
#' \describe{
#'   \item{list}{Common cluster between the three distance metrics}
#' }
#' @source source info (or `"manually created"`)
#' @examples
#' data(ex.common.clusters)
#' head(ex.common.clusters)
"ex.common.clusters"



#' Nodenames Example
#'
#' List of gene names in dataset.
#'
#' @format Dataframe of 10 rows and 1 column:
#' \describe{
#'   \item{Gene.Names}{Name of the gene}
#' }
#' @source "created by MakeDBInput"
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
#'   \item{Gene.1}{The gene interacting with Gene.2}
#'   \item{Gene.2}{The gene interacting with Gene.1}
#'   \item{Weight}{Strength of the interaction}
#'   \item{edgeType}{Type of interaction between Gene.1 and Gene.2}
#' }
#' @source "created by FindPPIEdges"
#' @examples
#' data(ex.ppi.network)
#' head(ex.ppi.network)
"ex.ppi.network"



#' Correlation Matrix Example
#'
#' Matrix relating the correlation values between ptms based on Spearman Similarity.
#'
#' @format Square matrix of 100 rows and columns:
#' \describe{
#'   \item{rows and columns}{Ptm names}
#'   \item{values}{Strength of correlation}
#' }
#' @source "created by MakeClusterList"
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
#'   \item{rows}{Post-translational modifications}
#'   \item{columns}{Experimental condition}
#'   \item{values}{Expression as found by mass spectrometry}
#' }
#' @source "manually created"
#' @examples
#' data(ex.ptmtable)
#' head(ex.ptmtable)
"ex.ptmtable"



#' TSNE matrices Example
#'
#' List of of lists of dataframes of ptm positions in 3D space to show their relationships with each other based on distance.
#'
#' @format List of 3 Distance metrics (Euclidean, Spearman, and the average of Spearman and Euclidean (SED)) and the clusters that form when they are run through tsne:
#' \describe{
#'   \item{Euclidean}{A list of 6 clusters of ptms based on their Euclidean distance and positional distance determined by tsne}
#'   \item{Spearman}{A list of 6 clusters of ptms based on their Spearman dissimilarity and positional distance determined by tsne}
#'   \item{SED}{A list of 6 clusters of ptms based on the average of their Euclidean distance and Spearman Dissimilarity and positional distance determined by tsne}
#' }
#' @source "created by MakeClusterList"
#' @examples
#' data(ex.clusters.list)
#' head(ex.clusters.list)
"ex.clusters.list"

