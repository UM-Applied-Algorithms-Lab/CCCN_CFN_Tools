##WORK IN PROG. JUST COPIED OVER FROM KEYFUNCTIONS##

#called :)
# Function to extract gene names from peptide names
pepgene <- function(peps) {
  unique(sapply(peps, function(x) unlist(strsplit(x, " ", fixed=TRUE))[1]))
}

#NOT CALLED
#' Extract Gene Names from Peptide Edge File
#'
#' This function extracts unique gene names from a peptide edge file.
#'
#' @param peptide.edgefile A data frame containing peptide edge information.
#'
#' @return A vector of unique gene names.
#' @export
#'
#' @examples
#' extract.gene.names(peptide.edgefile)
# Function to extract gene names from peptide edge file
extract.gene.names <- function(peptide.edgefile) {
  peps <- c(peptide.edgefile[,1], peptide.edgefile[,2])
  genes <- unique(sapply(peps, function(x) unlist(strsplit(x, " ", fixed=TRUE))[1]))
  return(genes)
}

#' Create Gene-Peptide Edges
#'
#' This function creates peptide edges for a given node list.
#'
#' @param nodelist A vector of node names.
#' @param pepkey A data frame containing peptide keys.
#'
#' @return A data frame of peptide edges with weights and edge types.
#' @export
#'
#' @examples
#' genepep.edges.3(nodelist, pepkey)
genepep.edges.3 <- function(nodelist, pepkey=ld.key) { #NEVER USED
  nodelist <- unique(nodelist)
  gpedges <- pepkey[pepkey$Gene.Name %in% nodelist, 1:2]
  names(gpedges)[1:2] <- c("Gene.1", "Gene.2")
  gpedges$edgeType <- "peptide"
  gpedges$Weight <- 1
  gpedges$Alt.Weight <- 100
  gpedges$Directed <- FALSE
  return(unique(gpedges))
}

#' Process Correlation Edges
#'
#' This function processes correlation edges from a given correlation matrix.
#'
#' @param cor_matrix A correlation matrix.
#' @param mode A string specifying the graph mode. Default is "lower".
#'
#' @return A data frame of correlation edges.
#' @export
#'
#' @examples
#' process_correlation_edges(cor_matrix)
# Function to process correlation edges
process_correlation_edges <- function(cor_matrix, mode="lower") {
  g <- igraph::graph_from_adjacency_matrix(as.matrix(cor_matrix), mode=mode, diag=FALSE, weighted="Weight")
  edges <- data.frame(igraph::as_edgelist(g))
  edges$Weight <- igraph::edge_attr(g)[[1]]
  edges$edgeType <- "correlation"
  edges$edgeType[edges$Weight <= -0.5] <- "negative correlation"
  edges$edgeType[edges$Weight >= 0.5] <- "positive correlation"
  edges <- edges[!is.na(edges$Weight),]
  names(edges)[1:2] <- c("Peptide.1", "Peptide.2")
  edges$Gene.1 <- sapply(edges$Peptide.1, pepgene)
  edges$Gene.2 <- sapply(edges$Peptide.2, pepgene)
  return(edges)
}

# Function to filter dual modifications
filter_dual_modifications <- function(edges, mod1, mod2) {
  dual_mod <- edges[intersect(grep(mod1, edges$Peptide.1), grep(mod2, edges$Peptide.2)), ]
  return(dual_mod)
}

# Function to analyze negative correlations
analyze_negative_correlations <- function(edges) {
  neg <- edges[edges$Weight < 0, ]
  vneg <- neg[abs(neg$Weight) >= 0.5, ] #Add commentMore actions #Syntax error?
  vvneg <- neg[abs(neg$Weight) > 0.543, ]

  neg_genes <- unique(neg$Gene.1)
  vneg_genes <- unique(vneg$Gene.1)
  vvneg_genes <- unique(vvneg$Gene.1)

  return(list(neg=neg, vneg=vneg, vvneg=vvneg,
              neg_genes=neg_genes, vneg_genes=vneg_genes, vvneg_genes=vvneg_genes))
}
