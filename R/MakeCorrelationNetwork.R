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

  #Triple loop to look through every combination of elements of the list and compare them
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
#' Make Correlation Network first finds the intersection between the Euclidean, Spearman, and SED cluster matrices in order to find
#' the intersection between the three groups. It then adds the Genes in these PTMs to a list of common clusters and turns it into an adjacency matrix.
#' This adjacency matrix is used to filter relevant data --- clusters --- from the Spearman correlation matrix. The resultant
#' cocluster correlation network shows strength of relationships between proteins using the common clusters between the three distance metrics.
#'
#' @param clusterlist A list of three-dimensional data frames used to represent ptms in space to show relationships between them based on distances. Based on Euclidean Distance, Spearman Dissimilarity, and SED (the average between the two)
#' @param ptm.correlation.matrix A data frame showing the correlation between ptms (as the rows and the columns). NAs are placed along the diagonal.
#' @param keeplength Only keep clusters of ptms whose size is larger than this parameter. (I.e keeplength = 2 then keep ("AARS", "ARMS", "AGRS") but not ("AARS", "ARMS"))
#' @param clusters.name Desired name for the common clusters output; defaults to common.clusters
#' @param cccn.name Desired name for the cocluster correlation network; defaults to cccn.matrix
#' @return The list of common clusters between all three distance metrics (Euclidean, Spearman, and SED) and a matrix showing strength of relationships between proteins using the common clusters between the three distance metrics (Euclidean, Spearman, and Combined (SED))
#' @export
#'
#' @examples
#' ex.ptm.cor <- ex.ptm.correlation.matrix
#' MakeCorrelationNetwork(ex.clusters.list, ex.ptm.cor, 1, "ex.clusters.common", "ex.cccn.matrix")
#' print(ex.clusters.common[c(1, 2, 3)])
#' utils::head(ex.cccn.matrix[, c(1,2,3,4,5)])
MakeCorrelationNetwork <- function(clusterlist, ptm.correlation.matrix, keeplength = 2, clusters.name = "clusters.common", filtered.ptm.cor.name = "filtered.cor.matrix", cccn.name = "cccn.matrix"){


  ### Helper fuction to take the submatrix from ptm.correlation.matrix of every row that starts with gene1 and every col that starts with gene2 ###
  correlation.value <- function(Gene1, Gene2){
    r <- filtered.ptm.correlation[
      grep(paste(Gene1, ""), rownames(filtered.ptm.correlation), value = TRUE), #Paste is required so that grep cannot find the gene in another gene. Such as, Gene1 = HAT will identify HIHATH as the same protein
      grep(paste(Gene2, ""), colnames(filtered.ptm.correlation), value = TRUE)] #Adding the space makes it so the entire word has to be there
    r <- as.matrix(r) #Needed? if singular value
    return(sum(r, na.rm = TRUE)) #Return sum
  }


  #### Generate the combined adjacency matrix by taking PTMs to Genes ###
  clusters.common <- FindCommonClusters(clusterlist[[1]], clusterlist[[2]], clusterlist[[3]], keeplength) #Call function at top of code
  genes.common <- c(clusters.common, recursive=TRUE) #Flatten the clusters
  filtered.ptm.cor <- ptm.correlation.matrix[genes.common, genes.common] #Filter the correlation Network by PTMs that cocluster

  gene.common <- sapply(genes.common, function(x) {unlist(strsplit(x, " ",  fixed=TRUE))[[1]]}) #Will just trim all elements for every subelement in a list of character vectors
  ulist <- unique(gene.common) #Use this for rownames and colnames

  cccn.matrix <- matrix(NA, nrow=length(ulist), ncol=length(ulist), dimnames=list(ulist, ulist)) #Initilize empty matrix

  for(d in 1:length(gene.common)){ #For every cluster
    cluster <- gene.common[[d]]    #Save the current cluster
    for(e in cluster){             #For every element in the cluster
      for(f in cluster){           #Connect E to F
        cccn.matrix[e, f] <- correlation.value(e, f) #This adds the correlation value
  }}}

  cccn.matrix[cccn.matrix==0] <- NA # Replace 0 with NA in the correlation matrix
  diag(cccn.matrix) <- NA # Remove self-loops by setting diagonal to NA


  ### Export Final Data Structure ###
  cccn.matrix[is.na(cccn.matrix)] <- 0 #Used to be function
  assign(clusters.name, clusters.common, envir = .GlobalEnv) #List of common clusters
  assign(cccn.name, cccn.matrix, envir = .GlobalEnv) #CoCluster Correlation Network
  assign(filtered.ptm.cor.name, filtered.ptm.cor, envir = .GlobalEnv) #Filtered correlation network


  ### Graphing ###
  graph <- igraph::graph_from_adjacency_matrix(cccn.matrix, mode = "lower", diag = FALSE, weighted = "Weight")
  plot(graph)
}
