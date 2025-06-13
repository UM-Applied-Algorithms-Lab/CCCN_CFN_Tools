#' Pathway Crosstalk Network
#' 
#' Converts Bioplanet pathways from (<https://tripod.nih.gov/bioplanet/>)  into a list of pathways whose elements are the genes in each pathway. Edge weights are either the PTM Cluster Weight or according to the Jaccard Similarity.
#' 
#' @param file Either the name of the bioplanet pathway .csv file OR the name of a dataframe loaded in environment, users should only pass in "yourfilename.csv"
#' @param clusterlist The list of coclusters made in MakeCorrelationNetwork
#' @return No clue (The result is a matrix with pathway names in columns and individual clusters as rows.)
#' @export
#'
#' @examples
#' print("TO DO")
PathwayCrosstalkNetwork <- function(file = "bioplanet.csv", clusterlist, PCNname = "pcn_matrix"){
#Read file in, converts to dataframe like with rows like: PATHWAY_ID | PATHWAY_NAME | GENE_ID | GENE_SYMBOL  
  #Loading .csv
  if(class(file) == "character"){
    if(!file.exists(file)) stop(paste(file, "not found. Plese check your working directory.")) #Error Catch
    bioplanet <- read.csv(file, stringsAsFactors = F) 
  }
  
  #Loading a dataframe for examples
  if(class(file) == "data.frame") bioplanet <- file
  
  #Turn bioplanet into a list of pathways. Pathways are character vectors comprised of gene names 
  pathways.list <- plyr::dlply(bioplanet, plyr::.(PATHWAY_NAME)) #Turn the bioplanet .csv into a list of data frames. Each data frame stores genes with the same PATHWAY_ID 
  pathways.list <- lapply(pathways.list, `[`, -c(1:3)) #Removes PATHWAY_ID | PATHWAY_NAME | GENE_ID from every single data frame in the list
  pathways.list <- lapply(pathways.list, unlist, use.names = FALSE) #Since data frames are 1 row, turn data frames into character vectors
  
  
#Jaccard Similarity#
  #Create a matrix whose [i, j] values are the intersection between the genes in pathway i and j
  matrix.jaccard <- matrix(0, nrow = length(pathways.list), ncol = length(pathways.list))
  
  #Rename
  rownames(matrix.jaccard) <- pathways.list
  colnames(matrix.jaccard) <- pathways.list

  #Populate matrix, diagonals must be 0 in order to prevent self-loops in 
  for (i in 1:length(pathways.list)) {
    if(i > length(pathways.list)) break #Or else out of bounds error will occur
    for (j in (i+1):length(pathways.list)) { #Populate upper triangular portion of the matrix
      p.intersect <- intersect(unlist(pathways.list[i]), unlist(pathways.list[j])) #Intersect
      p.union     <- union(unlist(pathways.list[i]), unlist(pathways.list[j]))     #Union 
      value <- length(p.intersect)/length(p.union) #Number of genes in intersect / Number of genes in union 
      if(value > 0) { #If value happens to be zero, just don't assign anything since NA is already there
        matrix.jaccard[i, j] <- value #Number of genes pathway i and j share
        matrix.jaccard[j, i] <- value #Since matrix is symmetrical
      }
    }}
  
  #Create igraph object using the matrix 
  pathways.graph <- igraph::graph_from_adjacency_matrix(matrix.jaccard, mode="lower", diag=FALSE, weighted="Weight")
    #plot(pathways.graph) #plot if desired
  
  #Interpret the jaccard matrix as an edgelist
  bioplanetjaccardedges <- data.frame(as_edgelist(pathways.graph)) #Convert to edgelist
  names(bioplanetjaccardedges) <- c("source", "target") #Rest is just renaming
  bioplanetjaccardedges$Weight <- edge_attr(pathways.graph)[[1]]
  bioplanetjaccardedges$interaction <- "pathway Jaccard similarity" 
  
  #BEFORE SAVING JACCARD SET DIAGION = 1
  
  
  
  
  
  
}