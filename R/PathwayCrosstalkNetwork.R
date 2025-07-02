# Cluster Pathway Evidence
#
# Helper function for obtaining cluster pathway evidence given a Cluster and Pathway
#
# @param genes A vector of strings containing all the gene names
# @param weights A vector of numerics containing all the weights for genes. Should map 1:1 to the genes vector. This is for tracking ambiguous PTMs and how they modify weights. 
# @param length  The length of the cluster (precomputed)
# @param pathway A pathway of genes/proteins
# @param p.list  A lookup table of pathways to calculate the # of pathways a gene, k, is in (pathways.list created in PathwayCrosstalkNetwork)
#
# @return A float value (cluster pathway evidence) between a cluster and pathway
ClusterPathwayEvidence <- function(genes, weights, length, pathway, p.lookup){
  sigma <- 0 #This is where the sum is stored
  
  if(is.null(weights)) weights <- rep(1, length(genes)) #if weights is not provided, assume no ambiguous PTMs (this should never be triggered in current code)
  
  #Calculate CPE score and add it to sigma
  for(k in 1:length(pathway)){
    temp <- sum(weights[pathway[k] == genes]) #The amount of times Gene k appears in cluster (can appear less than 1 time if Gene is apart of an ambiguous set like AARS ubi k747; ABCB1 p n40)
    temp <- temp / (p.lookup[pathway[k]]*length) #Divide temp by the number of times Gene k appears in pathways in the pathway list. Multiply # of time gene appears in all pathways by the size of the cluster as apart of the large cluster penalty.
    sigma <- sigma + temp #Add this value to sigma. 
  }
  return(sigma) #Return
}



# PTP Evidence Edge
# 
# Checks 2 columns of a matrix and returns the sum of all rows where both values in the row are nonzero
#
# @param pathway1 A column of CPE.matrix
# @param pathway2 A column of CPE.matrix
#
# @return The sum
PTP.evidence.edge <- function(col1, col2){
  sigma <- 0 #Zero, will return this
  for(i in 1:length(col1)){ #Both pathways will always be the same length
    if(col1[i] != 0 && col2[i] != 0) sigma <- sigma + col1[i] + col2[i] #If BOTH values in a row i from col1 and col2 are nonzero, add them both to sigma. If either one has a zero, add neither. 
  }
  return(sigma) #return
}



#' Pathway Crosstalk Network
#'
#' Converts Bioplanet pathways from (<https://tripod.nih.gov/bioplanet/>)  into a list of pathways whose elements are the genes in each pathway. Edge weights are either the PTM Cluster Weight or according to the Jaccard Similarity.
#'
#' @param file Either the name of the bioplanet pathway .csv file OR the name of a dataframe. Lines of bioplanet should possess 4 values in the order "PATHWAY_ID","PATHWAY_NAME","GENE_ID","GENE_SYMBOL". Users should only pass in "yourfilename.csv"
#' @param clusterlist The list of common clusters between all three distance metrics (Euclidean, Spearman, and SED) made in MakeCorrelationNetwork
#' @param edgelist.name Intended for use in Cytoscape. The desired name of the Pathway to Pathway edgelist file created ('.csv' will automatically be added to the end for you); defaults to edgelist
#' @return An edgelist file that is created in the working directory. Contains pathway source-target columns, with edge weights of their jaccard similarity and their Cluster-Pathway Evidence score
#' @export
#'
#' @examples
#' PathwayCrosstalkNetwork(ex.bioplanet, ex.common.clusters)
PathwayCrosstalkNetwork <- function(file = "bioplanet.csv", clusterlist, edgelist.name = "PTPedgelist"){
  #Read file in, converts to dataframe like with rows like: PATHWAY_ID | PATHWAY_NAME | GENE_ID | GENE_SYMBOL
  if(class(file) == "character"){   #Loading .csv
    if(!file.exists(file)) stop(paste(file, "not found. Plese check your working directory.")) #Error Catch
    bioplanet <- utils::read.csv(file, stringsAsFactors = F)
  }

  if(class(file) == "data.frame") bioplanet <- file   #Loading a dataframe for examples

  #Turn bioplanet into a list of pathways. Pathways are character vectors comprised of gene names
  pathways.list <- plyr::dlply(bioplanet, plyr::.(PATHWAY_NAME)) #Turn the bioplanet .csv into a list of data frames. Each data frame stores genes with the same PATHWAY_ID
  pathways.list <- lapply(pathways.list, `[`, 4) #Modifies all data frames in the list to only have the GENE_SYMBOL column. Uses [] as a function, which I did not know you could do in R. Very cool.
  pathways.list <- lapply(pathways.list, unlist, use.names = FALSE) #Since data frames are 1 row, turn data frames into character vectors



  ###Jaccard Similarity### - Pathways x pathways where value is the intersection divided by union
  jaccard.matrix <- matrix(0, nrow = length(pathways.list), ncol = length(pathways.list)) #Create an empty matrix that is pathways x pathways

  #Rename
  rownames(jaccard.matrix) <- names(pathways.list)
  colnames(jaccard.matrix) <- names(pathways.list)

  #Populate matrix using [i, j] values are the intersection/union between pathway i and j, diagonals must be 0 in order to prevent self-loops in
  for (i in 1:length(pathways.list)) {
    for (j in 1:length(pathways.list)) { #Populate symmetrical matrix
      p.intersect <- length(intersect(pathways.list[[i]], pathways.list[[j]])) #Length of Intersect
      p.union     <- length(pathways.list[[i]]) + length(pathways.list[[j]]) - p.intersect #Length of Union
      value <- p.intersect/p.union #Number of genes in intersect / Number of genes in union
      #If value happens to be zero, just don't assign anything since 0 is already there
      if(value > 0) jaccard.matrix[i, j] <- value #Number of genes pathway i and j share
    }}


  ###Pathway Cluster Evidence###  - A matrix of pathways x pathways whose values are found by using a custom formula that relates clusters and pathways

  CPE.matrix <- matrix(0, nrow = length(clusterlist), ncol = length(pathways.list))
  rownames(CPE.matrix) <- names(clusterlist) #Names
  colnames(CPE.matrix) <- names(pathways.list)

  #Populate Matrix
  pathways.lookup <- table(bioplanet$GENE_SYMBOL) #Create a lookup table for how many times each gene appears in the pathways list 
  
  
  for(a in 1:nrow(CPE.matrix)){
    #Precomputation steps for cluster
    #Use Gene names, NOT ptms. Note for anyone viewing these, data structures, names get messed up at this step due to R's c function
    gene.names <- sapply(clusterlist[[a]], function(x) strsplit(x, ";", fixed=TRUE)) #Turn all ambiguous proteins into a list which will be "Flattened out" in the next line
    gene.weights <- sapply(gene.names, function(y) rep(1/length(y), length(y) )) #Create weights for ambiguous PTMs that map onto gene.names,
    
    gene.names <- c(clusterlist[[a]], recursive=TRUE, use.names=FALSE) #Turns list of lists of lists into a character vector, easier to work with. Names will get messed up at this part
    gene.weights <- c(gene.weights, recursive=TRUE, use.names=FALSE) #Perform the same operation on weights to keep them mapped
    
    #Create a gene lookup table that stores the number of times every gene appears in the cluster. Accessed like: temp[['ABCA3']] gives the # of times ABCA3 appears. 
    gene.lookup <- table(sapply(gene.names, function(z) unlist(strsplit(z, " ",  fixed=TRUE))[1])) #Convert all PTMs to genes by cutting off modifications like "ubi 470"
    cluster.length <- length(clusterlist[[a]]) #Precompute the length
    
    
    #For every pathway for the given cluster, call ClusterPathwayEvidence (at top) for the CPE.matrix
    for(b in 1:ncol(CPE.matrix)) CPE.matrix[a, b] <- sum(gene.lookup[pathways.list[[b]]], na.rm=TRUE)/(sum(pathways.lookup[pathways.list[[b]]], na.rm=TRUE)*cluster.length)
  }



  ###Generate PCN network### -
  temp.rows <- apply(CPE.matrix, 1, function(x){colnames(CPE.matrix)[x!=0]}) #Mark valuable clusters in CPE.matrix. Creates a list of vectors that contain pathways connections where there is a nonzero weight. 1 Vector per row.
  if(length(temp.rows) == 0) stop("No Cluster Pathway Evidence found (Matrix is empty). Please ensure clusters and bioplanet have overlap.") #Error catch- Not worth continuing as a less helpful error will happen in the next line given a length of zero.
  temp.rows <- temp.rows[sapply(temp.rows, function(y){length(y)>=2})] #Remove every vector from temp.rows that below the length threshold (2)
  if(length(temp.rows) == 0) stop("No Cluster Pathway Evidence greater than 2 found") #Error catch- Not worth continuing as a less helpful error will happen in the next line given a length of zero.
  
  #Create data frame: Pathway to Pathway edgelist
  size <- sum(sapply(temp.rows, function(x) (length(x) * (length(x)-1))/2)) #This may look bad but it's just number of possible permutations where order doesnt matter bc I didn't want to import a package.
  PTP.edgelist <- data.frame(source = rep("-", size), target = rep("-", size), Jaccard_weight = rep(0, size), PTP_inter_evidence= rep(0, size))

  #Populate data frame
  track <- 1 #First Empty row in the data frame
  for(i in 1:length(temp.rows)){
    nodes <- combn(temp.rows[[i]], 2) #Get every node pair (permutations where order doesn't matter of a string vector). Stored as a matrix.
    cluster <- names(temp.rows)[[i]]  #Get the name of the cluster that the connection was found in
    for(pathway in asplit(nodes, 2)) { #Add all node pairings to data frame. This code splits the matrix that stores the permutations
      PTP.edgelist[track, 1:2] <- pathway #Add row from nodes to empty spot in the edgefiles
      path1 <- pathway[[1]] #These are strings representing the names of the pathways. Used to access data from data structures.
      path2 <- pathway[[2]]
      PTP.edgelist[track, 3] <- jaccard.matrix[path1, path2] #Add the jaccard weight to the edgelist
      PTP.edgelist[track, 4] <- PTP.evidence.edge(CPE.matrix[,path1], CPE.matrix[,path2]) #Pass in two columns of CPE.matrix to custom function. NOT SURE IF CORRECT. Take the sum of two CPE's, between the same cluster for pathway j[[1]] and j[[2]]

      track <- track+1 #Increase tracker
    }}


  ###Debug Variable Names### - DELETE ME
  assign("jaccard.matrix", jaccard.matrix, envir = .GlobalEnv) #DEBUG
  assign("CPE.matrix", CPE.matrix, envir = .GlobalEnv)         #DEBUG
  assign(edgelist.name, PTP.edgelist, envir = .GlobalEnv)      #DEBUG
  assign("temp.rows", temp.rows, envir = .GlobalEnv)           #DEBUG



  ###Save edgefile for cytoscape plotting###
  filename <- paste(edgelist.name, ".csv", sep="") #Name of the file created with .csv appended
  write.csv(PTP.edgelist, file = filename, row.names = FALSE) #Save to files for cytoscape... Correct formatting?

  #Tell the user where their files got put
  cat(filename, "made in directory:", getwd())
}
