# Cluster Pathway Evidence
#
# Helper function for obtaining cluster pathway evidence given a Cluster and Pathway
#
# @param cluster A cluster of genes/proteins (cannot be PTMs)
# @param pathway A pathway of genes/proteins
# @param p.list  A list of pathways to calculate the # of pathways a gene, k, is in (pathways.list created in PathwayCrosstalkNetwork)
#
# @return A float value representing cluster pathway evidence between a cluster and pathway
ClusterPathwayEvidence <- function(cluster, pathway, p.list){
  sigma <- rep(0, length(pathway)) #Cluster Pathway Evidence will be found by taking the sum of this vector

  #Use Gene names, NOT ptms. Note for anyone viewing these, data structures, names get messed up at this step due to R's c function
  cluster.format <- sapply(cluster, function(x) strsplit(x, "; ", fixed=TRUE)) #Turn all ambiguous proteins into a list which will be "Flattened out" in the next line
  cluster.weights <- sapply(cluster.format, function(y) rep(1/length(y), length(y) )) #Create weights for ambiguous PTMs that map onto cluster.format, 
  
  cluster.format <- c(cluster, recursive=TRUE, use.names=FALSE) #Turns list of lists of lists into a character vector, easier to work with. Names will get messed up at this part
  cluster.weights <- c(cluster.weights, recursive=TRUE, use.names=FALSE) #Perform the same operation on weights to keep them mapped
  cluster.format <- unique(sapply(cluster.format, function (z) unlist(strsplit(z, " ",  fixed=TRUE))[1])) #Convert all PTMs to genes by cutting off modifications like "ubi 470" and remove duplicates!

  #Calculate CPE score and add it to sigma
  for(k in 1:length(pathway)){
    temp <- sum(cluster.weights[pathway[k] == cluster.format]) #Numerator: The amount of times Gene k appears in cluster
    temp <- temp / sum(sapply(p.list, function(x) pathway[[k]] %in% x)) #Divide temp by the number of times Gene k appears in pathways in the pathway list
    sigma[k] <- temp #Assign this value to sigma[k]
  }

  #Large cluster penalty -> as cluster size increases, CPE decreases. This is applied to every element in sigma
  sigma <- sigma*(1/length(cluster))

  return(sum(sigma)) #Return
}


#' Pathway Crosstalk Network
#'
#' Converts Bioplanet pathways from (<https://tripod.nih.gov/bioplanet/>)  into a list of pathways whose elements are the genes in each pathway. Edge weights are either the PTM Cluster Weight or according to the Jaccard Similarity.
#'
#' @param file Either the name of the bioplanet pathway .csv file OR the name of a dataframe. Lines of bioplanet should possess 4 values in the order "PATHWAY_ID","PATHWAY_NAME","GENE_ID","GENE_SYMBOL". Users should only pass in "yourfilename.csv"
#' @param clusterlist The list of coclusters made in MakeCorrelationNetwork
#' @param edgelist The desired name for the edgelist file that PCN will make. Should NOT contain any file extensions like .csv, this step will add that for you. Intended for graphing in Cytoscape.
#' @export
#'
#' @examples
#' PathwayCrosstalkNetwork(ex.bioplanet, ex.common.clusters)
PathwayCrosstalkNetwork <- function(file = "bioplanet.csv", clusterlist, edgelist.name = "edgelist"){
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
  rownames(matrix.jaccard) <- names(pathways.list)
  colnames(matrix.jaccard) <- names(pathways.list)

  #Populate matrix, diagonals must be 0 in order to prevent self-loops in
  for (i in 1:length(pathways.list)) {
    for (j in 1:length(pathways.list)) { #Populate symmetrical matrix
      p.intersect <- length(intersect(pathways.list[[i]], pathways.list[[j]])) #Length of Intersect
      p.union     <- length(pathways.list[[i]]) + length(pathways.list[[j]]) - p.intersect #Length of Union
      value <- p.intersect/p.union #Number of genes in intersect / Number of genes in union
      #If value happens to be zero, just don't assign anything since 0 is already there
      if(value > 0) matrix.jaccard[i, j] <- value #Number of genes pathway i and j share
    }}

  #Create igraph object using the matrix
  pathways.graph <- igraph::graph_from_adjacency_matrix(matrix.jaccard, mode="lower", diag=FALSE, weighted="Weight")
  #plot(pathways.graph) #plot if desired

  #Interpret the jaccard matrix as an edgelist
  PCN.jaccardedges <- data.frame(igraph::as_edgelist(pathways.graph)) #Convert to edgelist
  names(PCN.jaccardedges) <- c("source", "target") #Rest is just renaming
  PCN.jaccardedges$Weight <- igraph::edge_attr(pathways.graph)[[1]]
  PCN.jaccardedges$interaction <- "pathway Jaccard similarity"


  ###Jaccard Edges assignment must be here in case of error throw in CPE step
  assign("PCN.jaccardedges", PCN.jaccardedges, envir = .GlobalEnv)

  ###Innit - Weights for Non-Ambiguous & Ambiguous PTMs###
  #gene.weights <- rep(1, length(do.call(c, clusterlist[[1]])) - length(do.call(c, clusterlist[[2]]))) #The weights of ALL non-ambiguous PTMs are 1. Number of non-ambiguous PTMs are found via # of Total ptms - # of ambiguous ptms
  #temp.ambig.weights <- unlist(sapply(clusterlist[[2]], function(x){return(rep(1/length(x), length(x)))})) #Weights of ambiguous PTMs are 1/# of ptms in the series. so Aars ubi k454; Abui Or an; Aars ubi k983 splits into 3 genes w/ weight of 1/3rd
  #gene.weights <- c(gene.weights, temp.ambig.weights) #Should look like a bunch of ones at the start followed by various weights < 1

  ###Generating pathway cluster evidence matrix###
  #A more correct way of doing things
  #MCN.data <- do.call(c, clusterlist[[1]], recursive = TRUE) #Clusters no longer matter - So convert into a list containing genenames. Called genevec in Mark's .rmd
  #geneweights <- rep(1, length(MCN.data)) #Create a vector of gene weights
  #genesinpathways <- MCN.data %in% pathways.genes #Which indices of geneweights cannot be 1.

  #My attempt at coding CPE formula - Will represent as a matrix; clusters x pathways
  CPE.Matrix <- matrix(0, nrow = length(clusterlist), ncol = length(pathways.list))
  rownames(CPE.Matrix) <- names(clusterlist) #Names
  colnames(CPE.Matrix) <- names(pathways.list)

  #Populate Matrix - TODO do NOT make CPE it's own function
  for(a in 1:nrow(CPE.Matrix)){
    for(b in 1:ncol(CPE.Matrix)){ #Use ClusterPathwayEvidence function (found at top)
      CPE.Matrix[a, b] <- ClusterPathwayEvidence(clusterlist[[a]], pathways.list[[b]], pathways.list) #Call Cluster Pathway Evidence (above)
  }}

  assign("CPE.Matrix", CPE.Matrix, envir = .GlobalEnv) #DEBUG
  ###Generate PCN network###

  #Isolate rows from CPE.Matrix
  temp.rows <- apply(CPE.Matrix, 1, function(x){colnames(CPE.Matrix)[x!=0]}) #Creates a list of vectors that contain pathways connections where there is a nonzero weight. 1 Vector per row.

  if(length(temp.rows) == 0) stop("No Cluster Pathway Evidence found (Matrix is empty). Please ensure clusters.common and bioplanet have overlap.") #Error catch- Not worth continuing as a less helpful error will happen in the loop given a length of zero.
  temp.rows <- temp.rows[sapply(temp.rows, function(y){length(y)>=2})] #Remove every vector from temp.rows that below the length threshold (2)
  #Create pathway x pathway submatrix with CPE weights here? 
  
  #Create data frame: Pathway to Pathway edgelist
  size <- sum(sapply(temp.rows, function(x) (length(x) * (length(x)-1))/2)) #This may look bad but it's just permutation where order doesnt matter bc I didn't want to import a package.
  PTP.edgelist <- data.frame(source = rep("-", size), target = rep("-", size), Jaccard_weight = rep(0, size), CPEweight = rep(0, size))

  #Populate data frame
  track <- 1 #Empty location in the data frame
  for(i in temp.rows){
    nodes <- combn(i, 2) #Get every node pair (permutations where order doesn't matter of a string vector)
    for(j in asplit(nodes, 2)) { #Add all node pairings to data frame
      PTP.edgelist[track, 1:2] <- j #Add row from nodes to empty spot in the edgefiles
      PTP.edgelist[track, 3] <- matrix.jaccard[j[[1]], j[[2]]] #Add the jaccard weight to the edgelist
      PTP.edgelist[track, 4] <- 1 #HOW DO YOU GET THE CPE Weights from CPE MATRIX???

      track <- track+1 #Increase tracker
    }}

  
  ###Assign Variable Names###
  assign("PTP.edgelist", PTP.edgelist, envir = .GlobalEnv) #DEBUG

  #Save edgefiles for cytoscape plotting
  filename <- paste(edgelist.name, ".csv", sep="") #Name of the file created with .csv appended
  write.csv(PTP.edgelist, file = filename, row.names = FALSE) #Save to files for cytoscape... Correct formatting?

  #Tell the user where their files got put
  cat(filename, "made in directory:", getwd())
}
