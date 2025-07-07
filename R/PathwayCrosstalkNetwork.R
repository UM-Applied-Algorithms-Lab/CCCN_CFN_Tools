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
  CPE.matrix <- matrix(NA, nrow = length(clusterlist), ncol = length(pathways.list))
  rownames(CPE.matrix) <- names(clusterlist) #Names
  colnames(CPE.matrix) <- names(pathways.list)

  #Populate Matrix
  pathways.temp <- as.data.frame(table(bioplanet$GENE_SYMBOL)) #Create a lookup table for how many times each gene appears in the pathways list 
  pathways.hash <- pathways.temp$Freq #Create a hashtable from the table, very important for lowering runtime
  names(pathways.hash) <- pathways.temp$Var1 #Now any string vector of genes like pathways.hash[[c("AARS", "ABCA1")]] will return the frequency of how many times those genes appear in the pathway list *in constant time*. Sum() to return the total  

  for(a in 1:nrow(CPE.matrix)){
    #Precomputation steps for cluster
    #Use Gene names, NOT ptms. Note for anyone viewing these, data structures, names get messed up at this step due to R's c function
    gene.names <- sapply(clusterlist[[a]], function(x) strsplit(x, ";", fixed=TRUE)) #Turn all ambiguous proteins into a list which will be "Flattened out" in the next line
    gene.weights <- sapply(gene.names, function(y) rep(1/length(y), length(y) )) #Create weights for ambiguous PTMs that map onto gene.names,
    
    gene.names <- c(clusterlist[[a]], recursive=TRUE, use.names=FALSE) #Turns list of lists of lists into a character vector, easier to work with. Names will get messed up at this part
    gene.weights <- c(gene.weights, recursive=TRUE, use.names=FALSE) #Perform the same operation on weights to keep them mapped
    
    #Create a gene lookup table that stores the number of times every gene appears in the cluster. Accessed like: temp[['ABCA3']] gives the # of times ABCA3 appears. 
    gene.temp <- as.data.frame(table(sapply(gene.names, function(z) unlist(strsplit(z, " ",  fixed=TRUE))[1]))) #Convert all PTMs to genes by cutting off modifications like "ubi 470"
    gene.hash <- gene.temp$Freq        #Create a hashtable from the gene.temp table
    names(gene.hash) <- gene.temp$Var1 #Same strategy as pathways.hash
    
    cluster.length <- length(clusterlist[[a]]) #Precompute the length
    
    #For every pathway for the given cluster, call ClusterPathwayEvidence (at top) for the CPE.matrix
    for(b in 1:ncol(CPE.matrix)){
      num <- sum(gene.hash[pathways.list[[b]]], na.rm=TRUE)
      dem <- (sum(pathways.hash[pathways.list[[b]]], na.rm=TRUE)*cluster.length)
      CPE.matrix[[a, b]] <- num/dem
    }
  }


  ###Generate PCN network### -for 
  temp.rows <- apply(CPE.matrix, 1, function(x){colnames(CPE.matrix)[x!=0]}) #Mark valuable clusters in CPE.matrix. Creates a list of vectors that contain pathways connections where there is a nonzero weight. 1 Vector per row.
  if(length(temp.rows) == 0) stop("No Cluster Pathway Evidence found (Matrix is empty). Please ensure clusters and bioplanet have overlap.") #Error catch- Not worth continuing as a less helpful error will happen in the next line given a length of zero.
  temp.rows <- temp.rows[sapply(temp.rows, function(y){length(y)>=2})] #Remove every vector from temp.rows that below the length threshold (2)
  if(length(temp.rows) == 0) stop("No Cluster Pathway Evidence greater than 2 found") #Error catch- Not worth continuing as a less helpful error will happen in the next line given a length of zero.
  
  
  

  ###Debug Variable Names### - DELETE ME
  assign("jaccard.matrix", jaccard.matrix, envir = .GlobalEnv) #DEBUG
  assign("CPE.matrix", CPE.matrix, envir = .GlobalEnv)         #DEBUG
  #assign(edgelist.name, PTP.edgelist, envir = .GlobalEnv)      #DEBUG
  assign("temp.rows", temp.rows, envir = .GlobalEnv)           #DEBUG



  ###Save edgefile for cytoscape plotting###
  #filename <- paste(edgelist.name, ".csv", sep="") #Name of the file created with .csv appended
  #utils::write.csv(PTP.edgelist, file = filename, row.names = FALSE) #Save to files for cytoscape... Correct formatting?

  #Tell the user where their files got put
  #cat(filename, "made in directory:", getwd())
}
