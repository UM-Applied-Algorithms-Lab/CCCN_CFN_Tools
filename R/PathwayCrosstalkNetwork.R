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
  if(is.character(file)){ #For loading .csv
    if(!file.exists(file)) stop(paste(file, "not found. Plese check your working directory.")) #Error Catch
    bioplanet <- utils::read.csv(file, stringsAsFactors = F) #Reads the file and turns it into a dataframe
  }

  if(is.data.frame(file)) bioplanet <- file   #For loading a dataframe for examples

  #Turn bioplanet into a list of pathways. Pathways are character vectors comprised of gene names
  pathways.list <- plyr::dlply(bioplanet, plyr::.(PATHWAY_NAME)) #Turn the bioplanet .csv into a list of data frames. Each data frame stores genes with the same PATHWAY_ID
  pathways.list <- lapply(pathways.list, `[`, 4) #Modifies all data frames in the list to only have the GENE_SYMBOL column. Uses [] as a function, which I did not know you could do in R. Very cool.
  pathways.list <- lapply(pathways.list, unlist, use.names = FALSE) #Since data frames are 1 row, turn data frames into character vectors





  ### Jaccard Similarity ### - List of edges, a data frame like PATHWAY | PATHWAY | JACCARD VALUE

  #Helper function to pass into an "apply" function for a matrix that contains character vectors. Matrix dimensions should be 2 columns x any number of rows
  find.jaccard.val <- function(charvectorrow){
    p.intersect <- length(intersect(charvectorrow[[1]], charvectorrow[[2]])) #Length of Intersect
    p.union     <- length(charvectorrow[[1]]) + length(charvectorrow[[2]]) - p.intersect #Length of Union
    value <- p.intersect/p.union#
    if(value == 0) return(NA)   #Return NA if 0 for igraph
    return(p.intersect/p.union) #Return the jaccard value
  }

  #Creating the edgelist
  PTPedgelist <- t(combn(names(pathways.list), 2)) #First two colums of the data frame like PATHWAY | PATHWAY; Found by taking the permutations of all pathway names in the string vectors
  combn.values <- combn(pathways.list, 2)  #Creates matrix whose rows should be passed into the above function
  jaccard.values <- apply(combn.values, 2, find.jaccard.val) #Find the jaccard value for each column in combn values
  PTPedgelist <- cbind(PTPedgelist, jaccard.values)#Attach the JACCARD VALUE column to PATHWAY | PATHWAY





  ###Pathway Cluster Evidence###  - A matrix of pathways x pathways whose values are found by using a custom formula that relates clusters and pathways
  CPE.matrix <- matrix(NA, nrow = length(clusterlist), ncol = length(pathways.list))
  rownames(CPE.matrix) <- names(clusterlist) #Names
  colnames(CPE.matrix) <- names(pathways.list)

  #Populate Matrix
  pathways.temp <- as.data.frame(table(bioplanet$GENE_SYMBOL)) #Create a lookup table for how many times each gene appears in the pathways list
  pathgene.count <- pathways.temp$Freq #Create a hashtable from the table, very important for lowering runtime
  names(pathgene.count) <- pathways.temp$Var1 #Now any string vector of genes like pathways.hash[[c("AARS", "ABCA1")]] will return the frequency of how many times those genes appear in the pathway list *in constant time*. Sum() to return the total

  for(a in 1:nrow(CPE.matrix)){
    #Precomputation steps for cluster
    #Use Gene names, NOT ptms. Note for anyone viewing these, data structures, names get messed up at this step due to R's c function
    gene.names <- sapply(clusterlist[[a]], function(x) strsplit(x, "; ", fixed=TRUE)) #Turn all ambiguous proteins into a list which will be "Flattened out" in the next line
    gene.weights <- sapply(gene.names, function(y) rep(1/length(y), length(y) )) #Create weights for ambiguous PTMs that map onto gene.names,

    gene.names <- c(gene.names, recursive=TRUE, use.names=FALSE) #Turns list of lists of lists into a character vector, easier to work with. Names will get messed up at this part
    gene.weights <- c(gene.weights, recursive=TRUE, use.names=FALSE) #Perform the same operation on weights to keep them mapped
    gene.names <- sapply(gene.names, function(z) unlist(strsplit(z, " ",  fixed=TRUE))[1]) #Crop off the PTM parts

    gene.count <- tapply(gene.weights, gene.names, sum) #Count the number of times genes appear in the vector, needs to be accessed FAST due to how many look ups is required.

    cluster.length <- length(clusterlist[[a]]) #Precompute the length

    #For every pathway for the given cluster, call ClusterPathwayEvidence (at top) for the CPE.matrix
    for(b in 1:ncol(CPE.matrix)){
      num <- gene.count[pathways.list[[b]]] #Calculate numerator - How many times each Protein from a pathway appears in the cluster (can appear multiple times due to PTMs, or less than 1 time(s) due to ambiguous PTMs)
      den <- pathgene.count[pathways.list[[b]]] #Calculate denominator - How many times each Protein in the pathway appears in the entire list of pathways * the length of the cluster
      value <- sum(num/den, na.rm=TRUE)*cluster.length #Vector wise division in num/den, then take the sum of the vector. Apply large cluster penalty
      if(value != 0) CPE.matrix[[a,b]] <- value #Important because I create a new vector by summing two columns. To see if two pathways have relation to the same cluster, I want to see if they have nonzero values in the same cluster (row). This is accomplished because only two nonzero values will result in an int, as int + NA = NA.
  }}






  ###Generate PCN network###
  CPE <- apply(PTPedgelist, 1, function(x) sum(CPE.matrix[,x[1]] + CPE.matrix[,x[2]], na.rm=TRUE)) #Get a vector of all the CPE weights for every permutation. For a CPE weight to be nonzero, they BOTH must have a mapping column in at least 1 area.
  CPE[CPE == 0] <- NA #Replace all 0s with NAs
  PTPedgelist <- cbind(PTPedgelist, CPE) #Bind all the columns together. Now Data structure is PATHWAY | PATHWAY | Jaccard | CPE
  PTPedgelist <- PTPedgelist[rowSums(is.na(PTPedgelist)) != 2, ] #Remove all rows that only have NA values for the jaccard and CPE values





  ###Debug Variable Names### - DELETE ME
  assign("CPE.matrix", CPE.matrix, envir = .GlobalEnv)   #DEBUG
  assign(edgelist.name, PTPedgelist, envir = .GlobalEnv) #DEBUG






  ###Save edgefile for cytoscape plotting###
  filename <- paste(edgelist.name, ".csv", sep="") #Name of the file created with .csv appended
  utils::write.csv(PTPedgelist, file = filename, row.names = FALSE) #Save to files for cytoscape... Correct formatting?

  #Tell the user where their files got put
  cat(filename, "made in directory:", getwd())
}
