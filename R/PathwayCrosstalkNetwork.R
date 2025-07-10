# PTPedgelist to igraph
#
# Helper function that converts a PTPedgelist into an igraph
#
# @param df Edgelist, a dataframe with first 2 columns as pathways x pathways, 3rd column as jaccard weights, 4th column as PTP weights. Meant to be used for the PTPedgelist created in Pathway Crosstalk Network
# @return An igraph object
PTPedgelist.to.igraph <- function(df){
  graph <- igraph::graph_from_data_frame(df[,1:2], directed = FALSE) #Create unweighted igraph
  igraph::E(graph)$weight <- df[,3] #Adds the jaccard value to the weights
  #igraph::E(graph)$weight <- ex.PTPedgelist[,4] #Adds the PTP value to weights (will overwrite jaccard)
  return(graph)
}

#' Pathway Crosstalk Network
#'
#' Converts Bioplanet pathways from (<https://tripod.nih.gov/bioplanet/>)  into a list of pathways whose elements are the genes in each pathway. Edge weights are either the PTM Cluster Weight or according to the Jaccard Similarity.
#'
#' @param file Either the name of the bioplanet pathway .csv file OR the name of a dataframe. Lines of bioplanet should possess 4 values in the order "PATHWAY_ID","PATHWAY_NAME","GENE_ID","GENE_SYMBOL". Users not well versed in R should only pass in "yourfilename.csv"
#' @param clusterlist The list of common clusters between all three distance metrics (Euclidean, Spearman, and SED). Can be made in MakeCorrelationNetwork
#' @param edgelist.name The desired name of the Pathway to Pathway edgelist file created ('.csv' will automatically be added to the end for you); defaults to edgelist. Intended for use in Cytoscape.
#' @return An edgelist file that is created in the working directory. Contains pathway source-target columns, with edge weights of their jaccard similarity and their Pathway-Pathway Evidence score
#' @export
#'
#' @examples
#' PathwayCrosstalkNetwork(ex.bioplanet, ex.common.clusters)
PathwayCrosstalkNetwork <- function(file = "bioplanet.csv", clusterlist, edgelist.name = "PTPedgelist"){


  #### Read file in, converts to dataframe like with rows like: PATHWAY_ID | PATHWAY_NAME | GENE_ID | GENE_SYMBOL ###
  if(is.character(file)){ #If Path to a .csv file (string input)
    if(!file.exists(file)) stop(paste(file, "not found. Plese check your working directory."))
    bioplanet <- utils::read.csv(file, stringsAsFactors = F) #Reads the file and turns it into a dataframe

  }else if(is.data.frame(file)){ bioplanet <- file #If data frame input

  }else stop(paste(class(file), "is not a supported file type. Please make sure you input a path to a .csv file or a data frame")) #If unsupported


  ### Turn bioplanet into a list of pathways. Pathways are character vectors comprised of gene names ###
  pathways.list <- plyr::dlply(bioplanet, plyr::.(PATHWAY_NAME)) #Turn the bioplanet .csv into a list of data frames. Each data frame stores genes with the same PATHWAY_ID
  pathways.list <- lapply(pathways.list, `[`, 4) #Modifies all data frames in the list to only have the GENE_SYMBOL column. Uses [] as a function, which I did not know you could do in R. Very cool.
  pathways.list <- lapply(pathways.list, unlist, use.names = FALSE) #Since data frames are 1 row, turn data frames into character vectors


  ### Functions to help caluclate the jaccard similarity ###
  find.jaccard.val <- function(charvectorrow){ #Helper function to pass into an "apply" function for a matrix that contains character vectors. Matrix dimensions should be 2 columns x any number of rows
    p.intersect <- length(intersect(charvectorrow[[1]], charvectorrow[[2]])) #Length of Intersect
    p.union     <- length(charvectorrow[[1]]) + length(charvectorrow[[2]]) - p.intersect #Length of Union
    jaccard.value <- p.intersect/p.union
    if(jaccard.value == 0) return(NA) #Return NA if 0 for igraph
    return(jaccard.value) #Return the jaccard value
  }

  combinations <- utils::combn #Mask the utils::combn with a more descriptive term


  ### Create the main data structure and add Jaccard Values ###
  PTPedgelist <- t(combinations(names(pathways.list), 2)) #Creating the edgelist. First two columns of the data frame like PATHWAY | PATHWAY; Found by finding all pairs pathway names
  combn.vector <- combinations(pathways.list, 2)  #The values of the above vector, performs the same operation to get pairs of string vectors containing gene names instead of pathway names. Prepares data for analysis
  jaccard.values <- apply(combn.vector, 2, find.jaccard.val) #Find the jaccard value for every pair of string vectors in combn.svector
  PTPedgelist <- cbind(PTPedgelist, jaccard.values) #Attach the JACCARD VALUE column to PATHWAY | PATHWAY

  assign("Jaccard.Full", PTPedgelist, envir = .GlobalEnv) #DEBUG - For viewing the full jaccard edgelist

  ### Pathway Cluster Evidence ###
  CPE.matrix <- matrix(NA, nrow = length(clusterlist), ncol = length(pathways.list)) #Initilize empty data structure, Clusters x Pathways
  rownames(CPE.matrix) <- names(clusterlist) #Names
  colnames(CPE.matrix) <- names(pathways.list)

  pathways.temp <- as.data.frame(table(bioplanet$GENE_SYMBOL)) #Create table for how many times each gene appears in the pathways list. Needs to be converted into a named vector for efficent runtime.
  pathgene.count <- pathways.temp$Freq #Transform into a named vector
  names(pathgene.count) <- pathways.temp$Var1 #Now any string vector of genes like pathways.hash[[c("AARS", "ABCA1")]] will return the frequency of how many times those genes appear in the pathway list *in constant time*. Sum() to return the total

  for(a in 1:nrow(CPE.matrix)){ #Populate the matrix - Perform the following steps for every cluster

    gene.names <- sapply(clusterlist[[a]], function(x) strsplit(x, "; ", fixed=TRUE)) #Turn all ambiguous proteins into a list within the list which will be "Flattened out" in the next line
    gene.weights <- sapply(gene.names, function(y) rep(1/length(y), length(y) )) #Create weights for ambiguous PTMs that map onto gene.names,

    gene.names <- c(gene.names, recursive=TRUE, use.names=FALSE) #Turns list of lists of lists into a character vector, easier to work with. Names will get messed up at this part
    gene.weights <- c(gene.weights, recursive=TRUE, use.names=FALSE) #Perform the same operation on weights to keep them mapped
    gene.names <- sapply(gene.names, function(z) unlist(strsplit(z, " ",  fixed=TRUE))[1]) #Crop the PTM part ("AARS ubi 474" -> "AARS")

    gene.count <- tapply(gene.weights, gene.names, sum) #Count the number of times genes appear in the vector, needs to be accessed FAST due to how many look ups is required.

    cluster.length <- length(clusterlist[[a]]) #Precompute

    for(b in 1:ncol(CPE.matrix)){ #For every pathway for the given cluster, perform the following steps

      num <- gene.count[pathways.list[[b]]] #Calculate numerator vector - How many times every Protein from a pathway appears in the cluster (can appear multiple times due to PTMs, or less than 1 time(s) due to ambiguous PTMs)
      den <- pathgene.count[pathways.list[[b]]] #Calculate denominator - How many times every Protein in the pathway appears in the entire list of pathways
      value <- sum(num/den, na.rm=TRUE)/cluster.length #Element wise division. Then take the sum of the vector. Apply large cluster penalty afterwards.
      if(value != 0) CPE.matrix[[a,b]] <- value #Important because I create a new vector by summing two columns. To see if two pathways have relation to the same cluster, I want to see if they have nonzero values in the same cluster (row). This is accomplished because only two nonzero values will result in an int, as int + NA = NA.
  }}


  ### Generate PCN network ###
  CPE <- apply(PTPedgelist, 1, function(x) sum(CPE.matrix[rowSums(!is.na(ex.CPE.matrix)) == 2,x]) #Get a vector of all the PTP weights for every pair of pathways using the CPE weights to filter. For a PTP weight to be non-NA, the PTP weight will be the sum of all clusters both pathways have nonzero CPEs in.
  CPE[CPE == 0] <- NA #Turn all 0s created in above line into NAs

  PTPedgelist <- cbind(PTPedgelist, CPE) #Bind all the columns together. Now Data structure is PATHWAY | PATHWAY | Jaccard | CPE
  PTPedgelist <- PTPedgelist[rowSums(is.na(PTPedgelist)) != 2, ] #Remove all rows that only have NA values for the jaccard and CPE values

  assign("CPE.matrix", CPE.matrix, envir = .GlobalEnv)   #Save for user viewing
  assign(edgelist.name, PTPedgelist, envir = .GlobalEnv)


  ### Save edgefile for cytoscape plotting ###
  filename <- paste(edgelist.name, ".csv", sep="") #Name of the file created with .csv appended
  utils::write.csv(PTPedgelist, file = filename, row.names = FALSE) #Save to files for cytoscape... Correct formatting?

  cat(filename, "made in directory:", getwd()) #Tell the user where their files got put
}
