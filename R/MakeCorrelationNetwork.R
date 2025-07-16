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
#' @param ptm.cccn.name The PTM correlation matrix filtered by PTMs that cocluster
#' @param gene.cccn.name Desired name for the cocluster correlation network; defaults to gene.cccn
#' @return The list of common clusters between all three distance metrics (Euclidean, Spearman, and SED) and a matrix showing strength of relationships between proteins using the common clusters between the three distance metrics (Euclidean, Spearman, and Combined (SED))
#' @export
#'
#' @examples
#' ex.ptm.cor <- ex.ptm.correlation.matrix
#' MakeCorrelationNetwork(ex.clusters.list, ex.ptm.cor, 1, "ex.clusters.common", "ex.ptm.cccn", "ex.gene.cccn")
#' print(ex.clusters.common[c(1, 2, 3)])
#' utils::head(ex.gene.cccn[, c(1,2,3,4,5)])
MakeCorrelationNetwork <- function(clusterlist, ptm.correlation.matrix, keeplength = 2, clusters.name = "clusters.common", ptm.cccn.name = "ptm.cccn", gene.cccn.name = "gene.cccn"){


  ### Helper fuction to take the submatrix from ptm.correlation.matrix of every row that starts with gene1 and every col that starts with gene2 ###
  correlation.value <- function(Gene1, Gene2){
    r <- ptm.cccn[
      grep(paste(Gene1, ""), rownames(ptm.cccn), value = TRUE), #Paste is required so that grep cannot find the gene in another gene. Such as, Gene1 = HAT will identify HIHATH as the same protein
      grep(paste(Gene2, ""), colnames(ptm.cccn), value = TRUE)] #Adding the space makes it so the entire word has to be there
    r <- as.matrix(r) #Needed? if singular value
    return(sum(r, na.rm = TRUE)) #Return sum
  }


  #### Generate the combined adjacency matrix by taking PTMs to Genes ###
  clusters.common <- FindCommonClusters(clusterlist[[1]], clusterlist[[2]], clusterlist[[3]], keeplength) #Call function at top of code
  genes.common <- c(clusters.common, recursive=TRUE) #Flatten the clusters
  ptm.cccn <- ptm.correlation.matrix[genes.common, genes.common] #Filter the correlation Network by PTMs that cocluster

  gene.common <- sapply(genes.common, function(x) {unlist(strsplit(x, " ",  fixed=TRUE))[[1]]}) #Will just trim all elements for every subelement in a list of character vectors
  ulist <- unique(gene.common) #Use this for rownames and colnames

  gene.cccn <- matrix(NA, nrow=length(ulist), ncol=length(ulist), dimnames=list(ulist, ulist)) #Initilize empty matrix

  for(d in 1:length(gene.common)){ #For every cluster
    cluster <- gene.common[[d]]    #Save the current cluster
    for(e in cluster){             #For every element in the cluster
      for(f in cluster){           #Connect E to F
        gene.cccn[e, f] <- correlation.value(e, f) #This adds the correlation value
  }}}

  gene.cccn[gene.cccn==0] <- NA # Replace 0 with NA in the correlation matrix
  diag(gene.cccn) <- NA # Remove self-loops by setting diagonal to NA


  ### Export Final Data Structure ###
  gene.cccn[is.na(gene.cccn)] <- 0 #Used to be function
  assign(clusters.name, clusters.common, envir = .GlobalEnv) #List of common clusters
  assign(gene.cccn.name, gene.cccn, envir = .GlobalEnv) #CoCluster Correlation Network
  assign(ptm.cccn.name, ptm.cccn, envir = .GlobalEnv) #Filtered correlation network


  ### Graphing ###
  graph <- igraph::graph_from_adjacency_matrix(gene.cccn, mode = "lower", diag = FALSE, weighted = "Weight")
  plot(graph)
}
