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

#' Create Adjacency Matrix
#'
#' This function creates an adjacency matrix for a given list element.
#'
#' @param list.element A list of elements to construct the adjacency matrix.
#'
#' @return A square matrix where rows and columns correspond to the input list elements.
#' @export
#'
#' @examples
#' MakeAdjMatrix(c("A", "B", "C"))
MakeAdjMatrix <- function(list.element) {
  list.el.mat <- matrix(1, nrow = length(list.element), ncol = length(list.element))
  rownames(list.el.mat) <- list.element
  colnames(list.el.mat) <- list.element
  return(list.el.mat)
}

#' Replace NAs with 0
#'
#' This function replaces all NA's in a data frame with 0 values.
#'
#' @param df A data frame where NA's are to be replaced with 0s.
#'
#' @return A data frame with NAs replaced by 0.
#' @export
#'
#' @examples
#' NA.to.zero.func(data.frame(a = c(NA, 1), b = c(2, NA)))
NA.to.zero.func <- function(df) {
  df0 <- df
  df0[is.na(df0)] <- 0
  return(df0)
}

#' Bind Matrices
#'
#' This function orders the adjacency matrix, finds the intersect with the correlation matrix (cccn_matrix), and puts NAs on the diagonal.
#'
#' @param cluster_list A list of clusters to generate adjacency matrices.
#' @param correlation_matrix A correlation matrix to align with the adjacency matrix.
#'
#' @return The CCCN matrix, the intersect between the correlation matrix and the ordered adjacency matrix.
#' @export
#'
#' @examples
#' BindMatrices(cluster_list, correlation_matrix)
BindMatrices <- function(cluster_list, correlation_matrix) {
  # Generate the combined adjacency matrix
  adj_matrix <- plyr::rbind.fill.matrix(plyr::llply(cluster_list, MakeAdjMatrix))
  rownames(adj_matrix) <- colnames(adj_matrix)
  
  # Order the adjacency matrix by row and column names
  adj_matrix_ordered <- adj_matrix[order(rownames(adj_matrix)), order(colnames(adj_matrix))]
  
  # Align the correlation matrix with the ordered adjacency matrix
  matched_rows <- intersect(rownames(adj_matrix_ordered), rownames(correlation_matrix)) #returns null for sample data
  matched_cols <- intersect(colnames(adj_matrix_ordered), colnames(correlation_matrix)) #returns null for sample data
  cccn_matrix <- correlation_matrix[matched_rows, matched_cols]
  
  # Replace NA values in the correlation matrix
  na_indices <- which(is.na(adj_matrix_ordered), arr.ind = TRUE)
  cccn_matrix <- replace(cccn_matrix, na_indices, NA)
  
  # Remove self-loops by setting diagonal to NA
  if (any(!is.na(diag(cccn_matrix)))) {
    diag(cccn_matrix) <- NA
  }
  
  # Return the adjacency and CCCN matrices as a list
  # the adj_matrix appears to be unused in future computations. PERHAPS DELETE?
  return(cccn_matrix)
}

#' Generate Correlation Network
#'
#' This function takes the cccn_matrix, replaces the NA values with zeroes, and turns it into an igraph object.
#'
#' @param cccn_matrix The CCCN matrix, the intersect between the correlation matrix and the ordered adjacency matrix produced from `BindMatrices`.
#'
#' @return An igraph object representing the correlation network.
#' @export
#'
#' @examples
#' CorrelationNetwork(cccn_matrix)
CorrelationNetwork <- function(cccn_matrix) {
  
  # Make igraph object, replacing NA with 0
  cccn_matrix0 <- NA.to.zero.func(cccn_matrix)
  graph <- igraph::graph_from_adjacency_matrix(as.matrix(cccn_matrix0), mode = "lower", diag = FALSE, weighted = "Weight")
  
  # Return the graph object
  return(graph)
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
  #Find common clusters
  list.common <- FindCommonClusters(eu_ptms_list, sp_ptms_list, sed_ptms_list, 1)
  
  #Find CCCN matrix
  cccn.matrix <- BindMatrices(list.common, dissimilarity.ptmtable)
  
  #Return correlation network 
  return(CorrelationNetwork(cccn.matrix))
}

