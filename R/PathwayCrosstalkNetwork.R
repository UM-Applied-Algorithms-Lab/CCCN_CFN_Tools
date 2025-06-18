#' Cluster Pathway Evidence
#' 
#' Helper function for obtaining cluster pathway evidence given a Cluster and Pathway
#' 
#' @param cluster A cluster of genes/proteins (cannot be PTMs)
#' @param pathway A pathway of genes/proteins
#' @param p.list  A list of pathways to calculate the # of pathways a gene, k, is in 
#' 
#' @return A float value representing cluster pathway evidence between a cluster and pathway
ClusterPathwayEvidence <- function(cluster, pathway, p.list){
  # Cluster Pathway Evidence = ∑( 1 /  # of pathways Gene k is in * size of cluster)
  #Currently just works with the denominator (Treats numerator as 1) as current iteration of the code cannot easily access the # of PTMs a Gene has
  sigma <- rep(0, length(pathway)) #Cluster Pathway Evidence will be found by taking the sum of this vector
  
  #Denominator for CPE function
  for(k in 1:length(pathway)) sigma[k] <- sum(sapply(p.list, function(x) pathway[[k]] %in% x)) #Assigns values to sigma for each K of the # of pathways that Gene k is in
  #NOTE TO SELF: Code should definitely be changed as even if this is a suitable solution (it is not) cluster should work way differently (Only used to *length??)
  sigma <- sigma * length(cluster) #Multiply by the length of cluster - Sigma should now contain denominators
  sigma <- 1/sigma #Divide all values of sigma by 1 (Should be PTMs from gene(i) pathway in cluster) to complete
  
  return(sum(sigma)) #Return the CPE
}

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
    bioplanet <- utils::read.csv(file, stringsAsFactors = F) 
  }
  
  #Loading a dataframe for examples
  if(class(file) == "data.frame") bioplanet <- file
  
  #Turn bioplanet into a list of pathways. Pathways are character vectors comprised of gene names 
  pathways.list <- plyr::dlply(bioplanet, plyr::.(PATHWAY_NAME)) #Turn the bioplanet .csv into a list of data frames. Each data frame stores genes with the same PATHWAY_ID 
  pathways.list <- lapply(pathways.list, `[`, -c(1:3)) #Removes PATHWAY_ID | PATHWAY_NAME | GENE_ID from every single data frame in the list
  pathways.list <- lapply(pathways.list, unlist, use.names = FALSE) #Since data frames are 1 row, turn data frames into character vectors
  
  
  ###Jaccard Similarity###
  #Create a matrix whose [i, j] values are the intersection/union between the genes in pathway i and j
  matrix.jaccard <- matrix(0, nrow = length(pathways.list), ncol = length(pathways.list))
  
  #Rename
  rownames(matrix.jaccard) <- pathways.list
  colnames(matrix.jaccard) <- pathways.list

  #Populate matrix, diagonals must be 0 in order to prevent self-loops in 
  for (i in 1:length(pathways.list)) {
    if(i > length(pathways.list)) break #Or else out of bounds error will occur
    for (j in (i+1):length(pathways.list)) { #Populate symmetrical matrix
      p.intersect <- length(intersect(unlist(pathways.list[i]), unlist(pathways.list[j]))) #Intersect
      p.union     <- length(pathways.list[i]) + length(pathways.list[i]) - p.intersect     #Union
      value <- p.intersect/p.union #Number of genes in intersect / Number of genes in union 
      if(value > 0) { #If value happens to be zero, just don't assign anything since 0 is already there
        matrix.jaccard[i, j] <- value #Number of genes pathway i and j share
        matrix.jaccard[j, i] <- value #Since matrix is symmetrical
      }
    }}
  
  #Create igraph object using the matrix 
  pathways.graph <- igraph::graph_from_adjacency_matrix(matrix.jaccard, mode="lower", diag=FALSE, weighted="Weight")
  #plot(pathways.graph) #plot if desired
  
  #Interpret the jaccard matrix as an edgelist
  jaccardedges <- data.frame(igraph::as_edgelist(pathways.graph)) #Convert to edgelist
  names(jaccardedges) <- c("source", "target") #Rest is just renaming
  jaccardedges$Weight <- igraph::edge_attr(pathways.graph)[[1]]
  jaccardedges$interaction <- "pathway Jaccard similarity" 
  
  
  ###Innit - Weights for Non-Ambiguous & Ambiguous PTMs###
  #gene.weights <- rep(1, length(do.call(c, clusterlist[[1]])) - length(do.call(c, clusterlist[[2]]))) #The weights of ALL non-ambiguous PTMs are 1. Number of non-ambiguous PTMs are found via # of Total ptms - # of ambiguous ptms 
  #temp.ambig.weights <- unlist(sapply(clusterlist[[2]], function(x){return(rep(1/length(x), length(x)))})) #Weights of ambiguous PTMs are 1/# of ptms in the series. so Aars ubi k454; Abui Or an; Aars ubi k983 splits into 3 genes w/ weight of 1/3rd
  #gene.weights <- c(gene.weights, temp.ambig.weights) #Should look like a bunch of ones at the start followed by various weights < 1
  
  
  ###Generating pathway cluster evidence###
  #A more correct way of doing things  
  #MCN.data <- do.call(c, clusterlist[[1]]) #Clusters no longer matter - So convert into a list containing genenames. Called genevec in Mark's .rmd
  #geneweights <- rep(1, length(MCN.data)) #Create a vector of gene weights
  #genesinpathways <- MCN.data %in% pathways.genes #Which indices of geneweights cannot be 1. 
  
  #My attempt at coding CPE formula - Will represent as a matrix; clusters x pathways
  CPE.Matrix <- matrix(0.000, nrow = length(list.found[[1]]), ncol = length(pathways.list))
  rownames(CPE.Matrix) <- list.found[[1]] #Naming, will look very ugly. 
  colnames(CPE.Matrix) <- pathways.list
  
  #Populate Matrix
  for(a in 1:nrow(CPE.Matrix)){
    for(b in 1:ncol(CPE.Matrix)){ #Use ClusterPathwayEvidence function (found at top)
      CPE.Matrix[a, b] <- ClusterPathwayEvidence(list.found$`Common Clusters`[[a]], pathways.list[[b]], pathways.list)
  }}
  
  View(CPE.Matrix)
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}