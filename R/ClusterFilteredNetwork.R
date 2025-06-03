##WORK IN PROG. JUST COPIED OVER FROM KEYFUNCTIONS##


#called :)
# Function to extract gene names from peptide names
pepgene <- function(peps) {
  unique(sapply(peps, function(x) unlist(strsplit(x, " ", fixed=TRUE))[1]))
}

#' Process Correlation Edges
#'
#' This function processes correlation edges from a given correlation matrix
#'
#' @param pBound The minimum value an edge weight can be to make it a positive correlation. Default = 0.5
#' @param nBound The maximum value an edge weight can be to make it a negative correlation. Default = -0.5
#' @param mode A string specifying the graph mode. Default is "lower".
#'
#' @return A data frame of correlation edges.
#' @export
#'
#' @examples
#' process_correlation_edges(cor_matrix)
# Function to process correlation edges
process_correlation_edges <- function(pBound = 0.5, nBound = -0.5, mode="lower") {
  #Formatting
  g <- igraph::graph_from_adjacency_matrix(cccn_matrix, mode=mode, diag=FALSE, weighted="Weight") #If igraph is needed, maybe have MCN "return" igraph instead?
  edges <- data.frame(igraph::as_edgelist(g)) #Turns igraph into edgelist into dataframe. (Couldn't we just use the correlation matrix? It should have the same info.) 
  edges$Weight <- igraph::edge_attr(g)[[1]]   #Add weights to edge data frame
  
  #Define positive or negative correlation
  edges$edgeType <- "correlation"
  edges$edgeType[edges$Weight >= pBound] <- "positive correlation" #If weight is greater/equal to pbound, mark correlation as positive
  edges$edgeType[edges$Weight <= nBound] <- "negative correlation" #If weight is less/equal to nbound, mark correlation as negative
  
  #
  edges <- edges[!is.na(edges$Weight),] #Corrects NA values, but...  
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
  vneg <- neg[abs(neg$Weight) >= 0.5, ]
  vvneg <- neg[abs(neg$Weight) > 0.543, ]
  
  neg_genes <- unique(neg$Gene.1)
  vneg_genes <- unique(vneg$Gene.1)
  vvneg_genes <- unique(vvneg$Gene.1)
  
  return(list(neg=neg, vneg=vneg, vvneg=vvneg,
              neg_genes=neg_genes, vneg_genes=vneg_genes, vvneg_genes=vvneg_genes))
}