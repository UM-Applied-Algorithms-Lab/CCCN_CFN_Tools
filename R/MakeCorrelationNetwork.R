# Helper function to find intersections of clusters
#
# Finds common elements between clusters in three lists.
#
# @param list1 A list of ptms.
# @param list2 A list of ptms.
# @param list2 A list of ptms.
# @param klength Minimum size of intersections to keep.
# @return A list containing, 1: A list of common clusters and 2: A character vector of removed, ambiguous PTMs .
FindCommonClusters <- function(list1, list2, list3, klength){

  #Convert lists into groups of ptms
  list1.ptms <- lapply(list1, function(x){x$"PTM.Name"}) #These are lists of character vectors
  list2.ptms <- lapply(list2, function(y){y$"PTM.Name"})
  list3.ptms <- lapply(list3, function(z){z$"PTM.Name"})

  returnme <- list() #Innit empty list

  #Triple loop to look through elements of the list and compare them
  for(a in 1:length(list1.ptms)){
    for(b in 1:length(list2.ptms)){
      for(c in 1:length(list3.ptms)){
        temp <- Reduce(intersect, list(list1.ptms[[a]], list2.ptms[[b]], list3.ptms[[c]])) #Take the intersection of 3 character vectors (as a vector)
        if(length(temp) > klength) returnme[[length(returnme)+1]] <- as.list(temp) #Add temp as a cluster (character vector) to returnme
  }}}
  #Return
  if(length(returnme) == 0) stop("No common clusters found") #This is for line 370, where the code will return out bounds error anyways if the list is empty!
  names(returnme) <- sapply(1:length(returnme), function(x){paste("Cluster", x)})
  return(returnme)
}

#' Make Correlation Network
#'
#' This function finds common clusters from the data derived in MakeClusterList and populates the correlation network (matrix) with values. The value of a cell in the correlation is the sum of a submatrix created from all PTMs of the row gene and all PTMs of the column gene.
#'
#' @param tsne.matrices List containing matrices that contain Euclidean, Spearman, and SED t-SNE coords respectively
#' @param ptm.correlation.matrix Correlation matrix made from ptm table
#' @param keeplength Only keep clusters of ptms whose size is larger than this parameter. (I.e keeplength = 2 then keep ("AARS", "ARMS", "AGRS") but not ("AARS", "ARMS"))
#' @param lists.name The desired name for the output of the list containing clusters of PTMs and Genes
#' @param cccn.name The desired name for the output of the Correlation Network Matrix
#' @export
#'
#' @examples
#' MakeCorrelationNetwork(ex.clusters.list, ex.ptm.correlation.matrix, 1, "ex.common.clusters", "ex.cccn")
MakeCorrelationNetwork <- function(clusterlist, ptm.correlation.matrix, keeplength = 2, clusters.name = "common.clusters", cccn.name = "cccn.matrix"){

  #Helper fuction to take the submatrix from ptm.correlation.matrix of every row that starts with gene1 and every col that starts with gene2
  correlation.value <- function(Gene1, Gene2){
    r <- ptm.correlation.matrix[
      grep(Gene1, rownames(ptm.correlation.matrix), value = TRUE),
      grep(Gene2, colnames(ptm.correlation.matrix), value = TRUE)]
    r <- as.matrix(r) #Unsure if needed
    return(sum(r, na.rm = TRUE)) #Return sum
  }

  #Find common clusters
  clusters.common <- FindCommonClusters(clusterlist[[1]], clusterlist[[2]], clusterlist[[3]], keeplength)

  # Generate the combined adjacency matrix by taking PTMs to Genes
  gene.common <- lapply(clusters.common, function(x) lapply(x,  function(y){unlist(strsplit(y, " ",  fixed=TRUE))[[1]]})) #Will just trim all elements for every subelement in a list of character vectors
  ulist <- unique(unlist(gene.common)) #Use this for rownames and colnames

  cccn.matrix <- matrix(NA, nrow=length(ulist), ncol=length(ulist), dimnames=list(ulist, ulist)) #Initilize empty matrix
  # Populate the empty matrix
  for(d in 1:length(gene.common)){ #For every cluster
    cluster <- gene.common[[d]]    #Save the current cluster
    for(e in cluster){             #For every element in the cluster
      for(f in cluster){           #Connect E to F
        cccn.matrix[e, f] <- correlation.value(e, f) #This adds the correlation value
  }}}

  # Replace 0 with NA in the correlation matrix
  cccn.matrix[cccn.matrix==0] <- NA

  # Remove self-loops by setting diagonal to NA
  if (any(!is.na(diag(cccn.matrix)))) diag(cccn.matrix) <- NA

  # Make igraph object, replacing NA with 0
  cccn.matrix[is.na(cccn.matrix)] <- 0 #Used to be function
  assign(clusters.name, clusters.common, envir = .GlobalEnv) #List of common clusters
  assign(cccn.name, cccn.matrix, envir = .GlobalEnv) #Matrix containing Euclidean t-SNE coords

  #Graphing
  graph <- igraph::graph_from_adjacency_matrix(cccn.matrix, mode = "lower", diag = FALSE, weighted = "Weight")
  plot(graph)
}
