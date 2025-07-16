#' Bioplanet Example
#'
#' A data frame version of a small subset of the Bioplanet.csv database file
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
#' head(ex.bioplanet)
"ex.bioplanet"



#' CoCluster Correlation Network
#'
#' Network of relationships between genes based on clustering of ptms from the three distance metrics.
#'
#' @keywords internal
#'
#' @format Square matrix with 19 rows and columns:
#' \describe{
#'   \item{rows/columns}{Gene}
#'   \item{values}{Strength of correlation}
#' }
#' @source "created by MakeCorrelationNetwork"
#' @examples
#' head(ex.cccn.matrix)
"ex.cccn.matrix"



#' Pathways List Example
#'
#' A named list with names of Pathways and string vectors that hold all the gene symbols from those pathways
#'
#' @keywords internal
#'
#' @format List of 1657 elements:
#' \describe{
#'   \item{pathway}{A pathway name with the list of genes involved}
#' }
#' @source "created from a Bioplanet.csv file in PathwayCrosstalkNetwork"
#' @examples
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
#' head(ex.gm.edges)
"ex.gm.edges"



#' Full GeneMANIA Network Example
#'
#' Network of relationships between genes based on information from the GeneMANIA database.
#'
#' @keywords internal
#'
#' @format Dataframe of 601,552 rows and 4 columns:
#' \describe{
#'   \item{Gene.1}{The gene interacting with Gene.2}
#'   \item{Gene.2}{The gene interacting with Gene.1}
#'   \item{raw.weights}{Strength of the interaction}
#'   \item{GM.weights}{Weights from Genemania file}
#' }
#' @source \url{https://apps.cytoscape.org/apps/genemania}
#' @examples
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
#'   \item{Ensembl.Gene.ID}{ID of the gene in Ensemble}
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
#' head(ex.gm.nodes)
"ex.gm.nodes"



#' List of Common Clusters Example
#'
#' Clusters that contain PTMs that are found to cocluster in all 3 T-SNE lists.
#'
#' @keywords internal
#'
#' @format List of 40 lists:
#' \describe{
#'   \item{list}{Common cluster between the three distance metrics}
#' }
#' @source "Created in MakeCorrelationNetwork, FindCommonClusters using 3 T-SNE lists"
#' @examples
#' head(ex.common.clusters)
"ex.common.clusters"



#' Nodenames Example
#'
#' List of gene names in dataset.
#'
#' @keywords internal
#'
#' @format Dataframe of 67 rows and 1 column:
#' \describe{
#'   \item{Gene.Names}{Name of the gene}
#' }
#' @source "created by MakeDBInput"
#' @examples
#' head(ex.nodenames)
"ex.nodenames"



#' PPI Network Example
#'
#' Network of relationships between genes based on PPI database STRINGdb.
#'
#' @keywords internal
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
#' head(ex.ppi.network)
"ex.ppi.network"



#' Correlation Matrix Example
#'
#' Matrix relating the correlation values between ptms based on Spearman Similarity.
#'
#' @keywords internal
#'
#' @format Square matrix of 497 rows and columns:
#' \describe{
#'   \item{rows and columns}{Ptm names}
#'   \item{values}{Strength of correlation}
#' }
#' @source "created by MakeClusterList"
#' @examples
#' head(ex.ptm.correlation.matrix)
"ex.ptm.correlation.matrix"



#' PTM Table Example
#'
#' Post-translational modification data under certain environmental conditions.
#'
#' @keywords internal
#'
#' @format Dataframe of 497 rows and 72 columns:
#' \describe{
#'   \item{first 3 columns}{Gene, modification, and Amino Acid respectively}
#'   \item{columns}{Names of drugs that were used on PTMs}
#'   \item{row}{Number of PTMs}
#'   \item{values}{Expression as found by mass spectrometry}
#' }
#' @source "manually created"
#' @examples
#' head(ex.ptmtable)
"ex.ptmtable"



#' Clusters List Example
#'
#' List of lists of dataframes containing PTMs that have been found to cluster based on their distance in 3D space for a specific distance metric
#'
#' @keywords internal
#'
#' @format List of 3 Distance metrics (Euclidean, Spearman, and the average of Spearman and Euclidean (SED)) and the clusters that form when they are run through tsne:
#' \describe{
#'   \item{Euclidean}{A list of 54 clusters of ptms based on their Euclidean distance and positional distance determined by tsne}
#'   \item{Spearman}{A list of 125 clusters of ptms based on their Spearman dissimilarity and positional distance determined by tsne}
#'   \item{SED}{A list of 39 clusters of ptms based on the average of their Euclidean distance and Spearman Dissimilarity and positional distance determined by tsne}
#' }
#' @source "created by MakeClusterList"
#' @examples
#' head(ex.clusters.list)
"ex.clusters.list"



#' Cluster Filtered Network Example
#'
#' Subset of the PPI network in the cocluster correlation network.
#'
#' @keywords internal
#'
#' @format Dataframe of 71 rows and 4 columns:
#' \describe{
#'   \item{Gene.1}{First gene interacting}
#'   \item{Gene.2}{Second gene interacting}
#'   \item{Interaction}{The type of interaction}
#'   \item{PPI.weight}{The weight of the interaction between these two genes as described by STRINGdb and GeneMANIA}
#' }
#' @source "created by MakeDBInput"
#' @examples
#' head(ex.cfn)
"ex.cfn"



#' STRINGdb Edges example
#'
#' List of edges between genes in dataset.
#'
#' @keywords internal
#'
#' @format Dataframe of 205 rows and 4 columns:
#' \describe{
#'   \item{Gene.1}{First gene interacting}
#'   \item{Gene.2}{Second gene interacting}
#'   \item{Interaction}{The type of interaction}
#'   \item{STRINGdb.combined_score}{The weight of the interaction between these two genes as described by STRINGdb}
#' }
#' @source "created by GetSTRINGdb"
#' @examples
#' head(ex.stringdb.edges)
"ex.stringdb.edges"



#' PTP Edgelist Example
#'
#' An edgelist filtered from the combination of every pathway based on nonzero Jaccard & PTP score
#'
#' @keywords internal
#'
#' @format Dataframe of 19 rows and 4 columns:
#' \describe{
#'   \item{V1}{A pathway}
#'   \item{V2}{A pathway}
#'   \item{jaccard.values}{The jaccard similarity of the two pathways sets of genes}
#'   \item{PTPscore}{The sum of all CPEs for every cluster where both pathways have a nonzero score}
#' }
#' @source "Created in PathwayCrosstalkNetwork"
#' @examples
#' head(ex.PTPedgelist)
"ex.PTPedgelist"

#' Funckey
#'
#' DESCRIPTION
#'
#' @keywords internal
#'
#' @format Dataframe of 24339 rows and 8 columns:
#' \describe{
#'   \item{Gene.Name}{Name of the gene}
#'   \item{Approved.Name}{Name of the protein}
#'   \item{Hugo.Gene.Family}{Family of the gene}
#'   \item{HPRD.Function}{Function of the protein}
#'   \item{nodeType}{Type of protein}
#'   \item{Domains}{Known domains of protein}
#'   \item{Compartment}{Location of the protein}
#'   \item{Compartment.Overview}{Alternate location of the protein}
#' }
#' @source ""
#' @examples
#' head(ex.funckey)
"ex.funckey"

#' CPE matrix example
#'
#' Scores that represents how likely it is that a pathway interacts with a given cluster. Matrix of Clusters x Pathways.
#'
#' @keywords internal
#' @format Matrix of 40 rows and 9 columns
#' \describe{
#'   \item{rows}{Clusters}
#'   \item{columns}{Pathways}
#'   \item{values}{Cluster Pathway Evidence}
#' }
#' @source "Made in PathwayCrosstalkNetwork"
#' @examples
#' ex.CPE.matrix[c(1,2) , c(1,2)]
"ex.CPE.matrix"

#' Filtered Correlation Matrix example
#'
#' A version of the Correlation Matrix filtered by coclustering PTMs
#'
#' @keywords internal
#' @format Matrix of 329 rows and columns
#' \describe{
#'   \item{rows}{Select PTMs}
#'   \item{columns}{Select PTMs}
#'   \item{values}{Same as in Correlation Matrix}
#' }
#' @source "Made in MakeCorrelationNetwork"
#' @examples
#' ex.ptm.cccn[c(1,2) , c(1,2)]
"ex.ptm.cccn"

#' All Tsne Coords
#'
#' A list of Tsne coords for Euclidian, Spearman, and SED
#'
#' @keywords internal
#' @format List of Euclidian, Spearman, and SED tsne coords
#' \describe{
#'   \item{First entry}{Data frame containing Euclidian Coords}
#'   \item{Second entry}{Data frame containing Spearman coords}
#'   \item{Third entry}{Data frame containing SED coords}
#' }
#' @source "Made in MakeClusterList"
#' @examples
#' ex.all.tsne.coords[1]
"ex.all.tsne.coords"
