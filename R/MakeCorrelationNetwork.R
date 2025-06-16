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

  #Find all the matching intersections of list1 and list2
  common <- c()    #Create an empty list to hold those intersections
  ambiguous <- c() #Hold onto ambiguous clusters, (ones with semicolons)
  a.names <- c()   #Stores the cluster an ambiguous ptm is taken from 
  c.names <- c()   #Stores the cluster and PTM a common cluster is taken from

  for(a in 1:length(list1.ptms)){ #Triple loop to look through elements of the list and compare them
    for(b in 1:length(list2.ptms)){
      for(c in 1:length(list3.ptms)){
        temp <- Reduce(intersect, list(list1.ptms[[a]], list2.ptms[[b]], list3.ptms[[c]])) #Take the intersection of 3 character vectors (as a vector)
        if(length(temp) > klength) {
          
          ### Handle Ambiguous Cluster ###
          ambiguous.PTMs <- grep("; ", temp) #Vector of indecies that contain ambigious PTMs
          a.temp <- temp[ambiguous.PTMs]
          temp.vec <- c()

          #Handling semicolons. Thank you Mark amazing code. I was struggling with this one
          for(i in a.temp){
            x <- unlist(strsplit(as.character(i), "; ", fixed=TRUE)) #Split semicolons
            genes <- (sapply(as.character(x),  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1])) #For new vectors split by semicolons, trim their names
            temp.vec <- c(temp.vec, unlist(genes)) #Add all the Genes to temp. Use a temp.vec if you want to seperate it by clusters
          }
          
          if(length(temp.vec) != 0) { #Only add to ambiguous if a ; is actually found
            ambiguous[[length(ambiguous)+1]] <- temp.vec #Add the temp vector to the ambiguous vectors list
            a.names <- c(a.names, paste("From Cluster", length(common)+1)) #Save name for what cluster the ambiguous ptm is taken from 
          }
          
          ### Handle Common Cluster ###
          if(length(ambiguous.PTMs) != 0) b.temp <- temp[-ambiguous.PTMs] #For non-ambiguous PTMs
          else b.temp <- temp #Really ugly solution but it works
          
          #Append solution
          names(b.temp) <- b.temp #Helps a bunch
          b.temp <- lapply(b.temp, function (x){unlist(strsplit(x, " ",  fixed=TRUE))[1]}) #Trim names from PTMs to Gene Names
          common[[length(common)+1]] <- b.temp  #And only add it to the list to return if it has enough values
  }}}}
  #Return
  if(length(common) == 0 && length(ambiguous) == 0) stop("No common clusters found") #This is for line 370, where the code will return out bounds error anyways if the list is empty!
  names(ambiguous) <- a.names #Rename (ugly solution but struggling to do it another way)
  returnme <- list(common, ambiguous)
  names(returnme) <- c("Common Genes", "Ambiguous PTMs")
  return(returnme)
}

#' Make Correlation Network
#' 
#' This function finds common clusters from the data derived in MakeClusterList and populates the correlation network (matrix) with values. The value of a cell in the correlation is the sum of a submatrix created from all PTMs of the row gene and all PTMs of the column gene.
#'
#' @param tsne.matrices #List containing matrices that contain Euclidean, Spearman, and SED t-SNE coords respectively
#' @param ptm.correlation.matrix #Correlation matrix made from ptm table
#' @param keeplength Only keep clusters of ptms whose size is larger than this parameter. (I.e keeplength = 2 then keep ["AARS", "ARMS", "AGRS"] but not ["AARS", "ARMS"])
#' @param lists.name #The desired name for the output of the list containing common and ambigious clusters
#' @param cccn.name #The desired name for the output of the Correlation Network Matrix
#' @export
#'
#' @examples
#' cccn.cfn.tools:::ex.MakeCorrelationNetwork(keeplength = 1)
MakeCorrelationNetwork <- function(tsne.matrices, ptm.correlation.matrix, keeplength = 2, lists.name = "list.found", cccn.name = "cccn_matrix"){

  #Helper fuction to take the submatrix from ptm.correlation.matrix of every row that starts with gene1 and every col that starts with gene2
  correlation.value <- function(Gene1, Gene2){
    r <- ptm.correlation.matrix[
      grep(Gene1, rownames(ptm.correlation.matrix), value = TRUE),
      grep(Gene2, colnames(ptm.correlation.matrix), value = TRUE)]
    r <- as.matrix(r) #Unsure if needed
    return(sum(r, na.rm = TRUE)) #Return sum
  }

  #Find common clusters
  list.found <- FindCommonClusters(tsne.matrices[[1]], tsne.matrices[[2]], tsne.matrices[[3]], keeplength)
  list.common <- list.found[[1]]

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

  # Replace 0 with NA in the correlation matrix
  cccn_matrix[cccn_matrix==0] <- NA

  # Remove self-loops by setting diagonal to NA
  if (any(!is.na(diag(cccn_matrix)))) diag(cccn_matrix) <- NA

  # Make igraph object, replacing NA with 0
  cccn_matrix[is.na(cccn_matrix)] <- 0 #Used to be function
  assign(lists.name, list.found, envir = .GlobalEnv) #List of common clusters
  assign(cccn.name, cccn_matrix, envir = .GlobalEnv) #Matrix containing Euclidean t-SNE coords

  #Graphing
  graph <- igraph::graph_from_adjacency_matrix(cccn_matrix, mode = "lower", diag = FALSE, weighted = "Weight")
  plot(graph)
}
