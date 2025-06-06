# Helper function to find intersections of clusters
#
# Finds common elements between clusters in three lists.
#
# @param list1 A list of ptms.
# @param list2 A list of ptms.
# @param list2 A list of ptms.
# @param klength Minimum size of intersections to keep.
# @return A list of common clusters.
FindCommonClusters <- function(list1, list2, list3, klength){

  #Convert lists into groups of ptms
  list1.ptms <- lapply(list1, function(x){x$"PTM.Name"}) #These are lists of character vectors
  list2.ptms <- lapply(list2, function(y){y$"PTM.Name"})
  list3.ptms <- lapply(list3, function(z){z$"PTM.Name"})

  #Find all the matching intersections of list1 and list2
  returnme <- list()  #Create an empty list to hold those intersections

  for(a in 1:length(list1.ptms)){ #Triple loop to look through elements of the list and compare them
    for(b in 1:length(list2.ptms)){
      for(c in 1:length(list3.ptms)){
        temp <- Reduce(intersect, list(list1.ptms[[a]], list2.ptms[[b]], list3.ptms[[c]])) #Take the intersection of 3 character vectors (as a vector)
        if(length(temp) > klength) {
          temp <- lapply(temp, function (x){ #Trim the names
            unlist(strsplit(x, " ",  fixed=TRUE))[1]})
          returnme[[length(returnme)+1]] <- unique(temp)  #And only add it to the list to return if it has enough values
  }}}}
  
  if(length(returnme) == 0) stop("No common clusters found") #This is for line 370, where the code will return out bounds error anyways if the list is empty!
  return(returnme)
}

#' Make Correlation Network
#'
#' This function creates a square adjacency matrix and aligns it with the correlation network and common clusters found from the data derrived in MakeClusterList,
#'
#' @param keeplength Only keep clusters of ptms whose size is larger than this parameter. (I.e keeplength = 2 then keep ["AARS", "ARMS", "AGRS"] but not ["AARS", "ARMS"])
#' @export
#'
#' @examples
#' MakeCorrelationNetwork(keeplength = 1)
MakeCorrelationNetwork <- function(keeplength = 2){

  #Helper fuction to take the submatrix from ptm.correlation.matrix of every row that starts with gene1 and every col that starts with gene2
  correlation.value <- function(Gene1, Gene2){
    r <- ptm.correlation.matrix[
      grep(Gene1, rownames(ptm.correlation.matrix), value = TRUE), 
      grep(Gene2, colnames(ptm.correlation.matrix), value = TRUE)]
    r <- as.matrix(r) #Unsure if neeeded
    return(sum(r, na.rm = TRUE)) #Return average
  }
  
  #Find common clusters
  list.common <- FindCommonClusters(eu_ptms_list, sp_ptms_list, sed_ptms_list, keeplength)

  # Generate the combined adjacency matrix
  ulist <- unique(unlist(list.common)) #Use this for rownames and colnames
  
  cccn_matrix <- matrix(NA, nrow=length(ulist), ncol=length(ulist), dimnames=list(ulist, ulist)) #Initilize empty matrix 
  # Populate the empty matrix
  for(d in 1:length(list.common)){ #For every cluster
    cluster <- list.common[[d]]    #Save the current cluster
    for(e in cluster){             #For every element in the cluster
      for(f in cluster){           #Connect E to F
        cccn_matrix[e, f] <- correlation.value(e, f) #This adds the correlation value 
  }}}
  
  ## FOR DEBUGGING ## 
  assign("DEBUG_matrix", cccn_matrix, envir = .GlobalEnv) 

  # Replace 0 with NA in the correlation matrix
  cccn_matrix[cccn_matrix==0] <- NA

  # Remove self-loops by setting diagonal to NA
  if (any(!is.na(diag(cccn_matrix)))) diag(cccn_matrix) <- NA

  # Make igraph object, replacing NA with 0
  cccn_matrix[is.na(cccn_matrix)] <- 0 #Used to be function
  assign("cccn_matrix", cccn_matrix, envir = .GlobalEnv) #Matrix containing Euclidean t-SNE coords
  
  graph <- igraph::graph_from_adjacency_matrix(cccn_matrix, mode = "lower", diag = FALSE, weighted = "Weight")
  plot(graph)
}