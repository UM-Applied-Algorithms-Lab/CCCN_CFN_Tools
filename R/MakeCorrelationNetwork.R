# Helper function to find intersections of clusters
#'
#' Finds common elements between clusters in three lists.
#'
#' @param list1 A list of ptms.
#' @param list2 A list of ptms.
#' @param list2 A list of ptms.
#' @param klength Minimum size of intersections to keep.
#' @return A list of common clusters.
#' @examples
#' FindCommonClusters(list1, list2, list3, keeplength)
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
        if(length(temp) > klength) returnme[[length(returnme)+1]] <- temp               #And only add it to the list to return if it has enough values
      }
    }
  }
  if(length(returnme) == 0) stop("No common clusters found") #This is for line 370, where the code will return out bounds error anyways if the list is empty!
  return(returnme)
}

#' Make Correlation Network
#' 
#' This functions runs all the steps required for make correlation network 
#' 
#' @param eu_ptms_list  Matrix containing Euclidean t-SNE coords
#' @param sp_ptms_list  Matrix containing Spearman t-SNE coords
#' @param sed_ptms_list Matrix containing combined t-SNE coords
#' 
#' @return An igraph object representing the correlation network.
#' @export
#' 
#' @examples
#' MakeCorrelationNetwork(eu_ptms_list, sp_ptms_list, sed_ptms_list, keeplength)
MakeCorrelationNetwork <- function(keeplength = 2){
  
  #This function creates a square adjacency matrix for a given list element.
  MakeAdjMatrix <- function(list.element) {
    list.el.mat <- matrix(1, nrow = length(list.element), ncol = length(list.element))
    rownames(list.el.mat) <- list.element
    colnames(list.el.mat) <- list.element
    return(list.el.mat)
  }
  
  #This function replaces all NA's in a data frame with 0 values.
  NA.to.zero.func <- function(df) {
    df0 <- df
    df0[is.na(df0)] <- 0
    return(df0)
  }
  
  #Find common clusters
  list.common <- FindCommonClusters(eu_ptms_list, sp_ptms_list, sed_ptms_list, 1)
  
  # Generate the combined adjacency matrix
  adj_matrix <- plyr::rbind.fill.matrix(plyr::llply(cluster_list, MakeAdjMatrix))
  rownames(adj_matrix) <- colnames(adj_matrix)
  
  #Rename the correlation matrix using ptmtable as a reference
  colnames(ptm.correlation.matrix) <- rownames(ptmtable)
  rownames(ptm.correlation.matrix) <- rownames(ptmtable)
  
  # Align the correlation matrix with the ordered adjacency matrix
  matched_rows <- intersect(rownames(adj_matrix), rownames(ptm.correlation.matrix)) #returns null for sample data
  matched_cols <- intersect(colnames(adj_matrix), colnames(ptm.correlation.matrix)) #returns null for sample data
  cccn_matrix  <- ptm.correlation.matrix[matched_rows, matched_cols]
  
  # Replace NA values in the correlation matrix
  na_indices <- which(is.na(adj_matrix), arr.ind = TRUE)
  cccn_matrix <- replace(cccn_matrix, na_indices, NA)
  
  # Remove self-loops by setting diagonal to NA
  if (any(!is.na(diag(cccn_matrix)))) diag(cccn_matrix) <- NA
  
  # Make igraph object, replacing NA with 0
  cccn_matrix0 <- NA.to.zero.func(cccn_matrix)
  Network <- igraph::graph_from_adjacency_matrix(as.matrix(cccn_matrix0), mode = "lower", diag = FALSE, weighted = "Weight")
  
  #Return correlation network 
  return(Network)
}

