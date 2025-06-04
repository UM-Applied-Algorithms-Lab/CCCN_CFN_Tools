##WORK IN PROG. JUST COPIED OVER FROM KEYFUNCTIONS##

#' ClusterFilteredNetwork
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
#' ClusterFilteredNetwork()
# Function to process correlation edges
ClusterFilteredNetwork <- function(pBound = 0.5, nBound = -0.5, mode="lower") {
  #Formatting
  g <- igraph::graph_from_adjacency_matrix(cccn_matrix, mode=mode, diag=FALSE, weighted="Weight") #If igraph is needed, maybe have MCN "return" igraph instead?
  edges <- data.frame(igraph::as_edgelist(g)) #Turns igraph into edgelist into dataframe. (Couldn't we just use the correlation matrix? It should have the same info.) 
  edges$Weight <- igraph::edge_attr(g)[[1]]   #Add weights to edge data frame
  
  #Define positive or negative correlation
  edges$edgeType <- "correlation"
  edges$edgeType[edges$Weight >= pBound] <- "positive correlation" #If weight is greater/equal to pbound, mark correlation as positive
  edges$edgeType[edges$Weight <= nBound] <- "negative correlation" #If weight is less/equal to nbound, mark correlation as negative
  
  #Rename
  names(edges)[1:2] <- c("Peptide.1", "Peptide.2") #Are these supposed to be specific names? 
  
  #Analyze negative correlations - DO THESE NEED TO BE USER DECIDED?? 
  neg <- edges[edges$Weight < 0, ]    #Get negative weights
  vneg <- neg[neg$Weight <= nBound, ] #Get negative correlative weights
  vvneg <- neg[neg$Weight < -0.543, ]  #Some other parameter? 
  
  #BROKEN - This implies genepep.edges.3 call. 
  neg_genes <- unique(neg$Gene.1)
  vneg_genes <- unique(vneg$Gene.1)
  vvneg_genes <- unique(vvneg$Gene.1)
  
  #Output
  print("Edges:")
  print(edges)
}
