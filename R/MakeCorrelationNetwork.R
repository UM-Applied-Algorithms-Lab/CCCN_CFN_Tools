#' Make Correlation Network
#'
#' Make Correlation Network filters a correlation matrix of PTMs by specific PTMs.
#' It groups the PTM correlation matrix based on the Genes of PTMs.
#' By summing these submatrices, it also produces a gene by gene cocluster correlation network shows strength of relationships between proteins using the common clusters between the three distance metrics.
#'
#' @param common.clusters A list of clusters. Ideally the ones found by MakeClusterList in common_clusters
#' @param ptm.correlation.matrix A data frame showing the correlation between ptms (as the rows and the columns). NAs are placed along the diagonal.
#' @param ptm.cccn.name Desired name of the PTM correlation matrix filtered by specific PTMs from input
#' @param gene.cccn.name Desired name for the cocluster correlation network; defaults to gene.cccn
#' @return The list of common clusters between all three distance metrics (Euclidean, Spearman, and SED) and a matrix showing strength of relationships between proteins using the common clusters between the three distance metrics (Euclidean, Spearman, and Combined (SED))
#' @export
#'
#' @examples
#' ex.ptm.cor <- ex.ptm.correlation.matrix
#' MakeCorrelationNetwork(ex.common.clusters, ex.ptm.cor, "ex.ptm.cccn", "ex.gene.cccn")
#' ex.ptm.cccn[1:5, 1:5]
#' ex.gene.cccn[1:5, 1:5]
MakeCorrelationNetwork <- function(common.clusters, ptm.correlation.matrix, ptm.cccn.name = "ptm.cccn", gene.cccn.name = "gene.cccn"){

  ### Helper fuction to take the submatrix from ptm.correlation.matrix of every row that starts with gene1 and every col that starts with gene2 ###
  correlation.value <- function(Gene1, Gene2){
    r <- ptm.cccn[
      grep(paste(Gene1, "."), rownames(ptm.cccn), value = TRUE), #Paste is required so that grep cannot find the gene in another gene. Such as, Gene1 = HAT will identify HIHATH as the same protein
      grep(paste(Gene2, ""), colnames(ptm.cccn), value = TRUE)] #Adding the space makes it so the entire word has to be there
    r <- as.matrix(r) #Needed? if singular value
    return(sum(r, na.rm = TRUE)) #Return sum
  }


  #### Generate the combined adjacency matrix by taking PTMs to Genes ###
  ptms.vector <- c(common.clusters, recursive=TRUE) #Flatten the clusters
  ptm.cccn <- ptm.correlation.matrix[ptms.vector, ptms.vector] #Filter the correlation Network by PTMs that cocluster

  gene.vector <- sapply(ptms.vector, function(x) {unlist(strsplit(x, " ",  fixed=TRUE))[[1]]}) #Will just trim all elements for every subelement in a list of character vectors
  ulist <- unique(gene.vector) #Use this for rownames and colnames

  gene.cccn <- matrix(NA, nrow=length(ulist), ncol=length(ulist), dimnames=list(ulist, ulist)) #Initilize empty matrix

  for(a in 1:length(ulist)){ #For every cluster
    for(b in 1:length(ulist)){
      gene.cccn[a, b] <- correlation.value(ulist[[a]], ulist[[b]]) #This adds the correlation value
    }
  }

  gene.cccn[gene.cccn==0] <- NA # Replace 0 with NA in the correlation matrix
  diag(gene.cccn) <- NA # Remove self-loops by setting diagonal to NA


  ### Export Final Data Structure ###
  gene.cccn[is.na(gene.cccn)] <- 0 #Used to be function
  assign(gene.cccn.name, gene.cccn, envir = .GlobalEnv) #CoCluster Correlation Network
  assign(ptm.cccn.name, ptm.cccn, envir = .GlobalEnv) #Filtered correlation network


  ### Graphing ###
  graph <- igraph::graph_from_adjacency_matrix(gene.cccn, mode = "lower", diag = FALSE, weighted = "Weight")
  plot(graph)
}
