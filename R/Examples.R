#' Make Cluster List Example
#' @keywords internal
ex.MakeClusterList <- function(ptmtable, toolong = 3.5){

  #SPEARMAN CALCULATION

  # Add if statement here to make sure functions are formatted correctly #
  # Ensure ptmtable is a data frame with numeric values #
  ptmtable.sp <- as.data.frame(lapply(ptmtable, as.numeric))

  # Calculate Spearman correlation #
  ptm.correlation.matrix <- stats::cor(t(ptmtable.sp), use = "pairwise.complete.obs", method = "spearman")

  # Replace diagonal with NA #
  diag(ptm.correlation.matrix) <- NA

  # Calculate dissimilarity #
  sp.diss.matrix <- 1 - abs(ptm.correlation.matrix)

  # Handle any remaining NA values by setting them to the maximum dissimilarity #
  max_dist_eu <- max(sp.diss.matrix, na.rm = TRUE)
  sp.diss.matrix[is.na(sp.diss.matrix)] <- max_dist_eu

  # Fix names of correlation matrix
  colnames(ptm.correlation.matrix) <- ptmtable$PTM
  rownames(ptm.correlation.matrix) <- ptmtable$PTM

  # Fix names of dissimilarity matrix
  colnames(sp.diss.matrix) <- ptmtable$PTM
  rownames(sp.diss.matrix) <- ptmtable$PTM

  # Run t-SNE #
  tsne_results <- Rtsne::Rtsne(sp.diss.matrix, dims = 3, perplexity = 15, theta = 0.25, max_iter = 5000, check_duplicates = FALSE, pca = FALSE)
  # Return t-SNE results #
  spearman_result = tsne_results$Y


  #EUCLIDEAN CALCULATION

  # Add if statement here to make sure functions are formatted correctly #
  # Convert the dataframe to a distance matrix using Euclidean distance #
  ptmtable.dist = as.matrix(stats::dist(ptmtable, method = "euclidean"))

  # Compute the maximum distance in the matrix, excluding NA values #
  max_dist = max(ptmtable.dist, na.rm = TRUE)

  # Replace NA values in the distance matrix with 100 times the maximum distance #
  ptmtable.dist[is.na(ptmtable.dist)] <- 100 * max_dist

  # Normalize the distance matrix by scaling it to a range from 0 to 100 #
  eu.diss.calc <- 100 * ptmtable.dist / max_dist

  # Apply t-SNE to the distance matrix to reduce dimensions to 3 #
  # Parameters: dims = 3 (3D output), perplexity = 15, theta = 0.25 (speed/accuracy trade-off) #
  # max_iter = 5000 (number of iterations), check_duplicates = FALSE (treat rows as unique) #
  # pca = FALSE (no initial PCA) #
  eu.ptms.tsne.list <- Rtsne::Rtsne(as.matrix(eu.diss.calc), dims = 3, perplexity = 15, theta = 0.25, max_iter = 5000, check_duplicates = FALSE, pca = FALSE)

  # Extract the t-SNE results from the output list #
  euclidean_result <- eu.ptms.tsne.list$Y

  #COMBINED CALCULATION

  #fix spearman thing; so do the exact same thing but no absolute value
  sp.diss.calc <- 1 - ptm.correlation.matrix                        # range goes from (-1 to 1) to (0 to 2)
  max_diss_sp <- max(sp.diss.calc, na.rm = TRUE)                    # find the max value (around 2)
  sp.diss.calc <- sp.diss.calc * (max_dist / max_diss_sp)           # SCALING. THIS IS WHERE SCALING OCCURS. All the values are scaled so biggest sp = biggest eu
  sp.diss.calc[is.na(sp.diss.calc)] <- 50 * max_dist_eu             # make the NAs roughly equal to 100
  sp.diss.calc <- as.matrix(sp.diss.calc)                           # turn into a matrix

  #fix euclidean rq
  eu.diss.calc <- as.matrix(eu.diss.calc)

  #find average
  combined_distance <- (sp.diss.calc + eu.diss.calc) / 2
  # Perform t-SNE on the combined distances #
  tsne_result <- Rtsne::Rtsne(as.matrix(combined_distance), dims = 3, perplexity = 15, theta = 0.25, check_duplicates = FALSE, pca = FALSE)
  sed_result <- tsne_result$Y


  #Nested function to analyze result
  clustercreate <- function(result, ptmtable){

    #Compute the minimum spanning tree connecting the points
    tsne.span2 <- vegan::spantree(stats::dist(result), toolong=toolong)

    #Find clusters that are connected based on toolong (distance?)
    result.disc2 <-  vegan::distconnected(stats::dist(result), toolong = toolong, trace = TRUE)  # test
    cat ("threshold dissimilarity", toolong, "\n", max(result.disc2), " groups","\n")

    #Create a plot of the clusters using vegan
    vegan::ordiplot(result)
    #lines(tsne.span2, result) #???
    vegan::ordihull(result, result.disc2, col="red", lwd=2)

    #Format a data frame
    result.span2.df <- data.frame(ptmtable$PTM)
    names(result.span2.df) <- "PTM.Name"
    result.span2.df$group <- result.disc2 #Add groups found above to the data frame

    #Convert data frame into a list of clusters (check doesn't like group but it's a column name)
    result.span2.list <- plyr::dlply(result.span2.df, plyr::.(group))  # GROUP LIST  !
    return(result.span2.list)

  } #END of nested function

  #Assign different analysises to global enviroment
  assign("ex.eu_ptms_list", clustercreate(euclidean_result, ptmtable), envir = .GlobalEnv) #Matrix containing Euclidean t-SNE coords
  assign("ex.sp_ptms_list", clustercreate(spearman_result, ptmtable), envir = .GlobalEnv)  #Matrix containing Spearman t-SNE coords
  assign("ex.sed_ptms_list", clustercreate(sed_result, ptmtable), envir = .GlobalEnv)      #Matrix containing combined t-SNE coords
  assign("ex.ptm.correlation.matrix", ptm.correlation.matrix, envir = .GlobalEnv)          #Correlation Matrix for later use
}





#' Make Correlation Network Example
#' @keywords internal
ex.MakeCorrelationNetwork <- function(keeplength = 2){

  #Helper fuction to take the submatrix from ptm.correlation.matrix of every row that starts with gene1 and every col that starts with gene2
  correlation.value <- function(Gene1, Gene2){
    r <- ex.ptm.correlation.matrix[
      grep(Gene1, rownames(ex.ptm.correlation.matrix), value = TRUE),
      grep(Gene2, colnames(ex.ptm.correlation.matrix), value = TRUE)]
    r <- as.matrix(r) #Unsure if neeeded
    return(sum(r, na.rm = TRUE)) #Return average
  }

  #Find common clusters
  list.common <- FindCommonClusters(ex.eu_ptms_list, ex.sp_ptms_list, ex.sed_ptms_list, keeplength)

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

  graph <- igraph::graph_from_adjacency_matrix(cccn_matrix, mode = "lower", diag = FALSE, weighted = "Weight")
  plot(graph)
}





########## FIX THIS ONE IDK WHAT IS REALLY HAPPENING WHEN THE EXAMPLE IS CALLED ##########

#' Make GM Input Example
#' @keywords internal
ex.make_db_input <- function(ex.cccn_matrix) {
  cat(ex.nodenames[[1]], sep = "\n")
}





#' Find PPI Edges Example
#' @keywords internal
ex.find_ppi_edges <- function(ex.cccn_matrix, db_filepaths = c()) {

  # Initialize the STRING database object
  string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=0, link_data="detailed", input_directory="")

  # Retrieve the proteins from the STRING database
  string_proteins <- string_db$get_proteins()
  print(dim(string_proteins))

  if (!"Gene.Names" %in% colnames(nodenames)) {
    stop("Column 'Gene.Names' not found in nodenames.")
  }

  # Map the genes to STRING IDs
  #please note that nodenames replaces the previous "input_dataset"; nodenames appears to work well :)
  #Gene.Names also replaces experimental
  mapped_genes <- string_db$map(nodenames, "Gene.Names", removeUnmappedRows = TRUE)
  print(utils::head(mapped_genes))

  # Retrieve the interactions for the mapped genes
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)

  # Convert protein IDs to gene names
  interactions$Gene.1 <- sapply(interactions$from, function(x) string_proteins[match(x, string_proteins$protein_external_id), "preferred_name"])
  interactions$Gene.2 <- sapply(interactions$to, function(x) string_proteins[match(x, string_proteins$protein_external_id), "preferred_name"])

  # Filter interactions based on evidence types
  #These were picked based on my limited knowledge of ppis
  #may easily be changed
  str.e <- interactions[interactions$experimental > 0, ]
  str.n <- interactions[interactions$neighborhood > 0, ]
  str.d <- interactions[interactions$database > 0, ]
  str.co <- interactions[interactions$cooccurence > 0, ]
  str.ce <- interactions[interactions$coexpression > 0, ]

  # Combine filtered interactions
  combined_interactions <- unique(rbind(str.e, str.n, str.d, str.co, str.ce))

  # Assign edge types
  combined_interactions$edgeType <- "STRINGdb"
  combined_interactions[combined_interactions$database > 0, "edgeType"] <- "database"
  combined_interactions[combined_interactions$neighborhood > 0, "edgeType"] <- "neighborhood"
  combined_interactions[combined_interactions$experimenal > 0, "edgeType"] <- "experimental"
  combined_interactions[combined_interactions$cooccurence > 0, "edgeType"] <- "cooccurence"
  combined_interactions[combined_interactions$coexpression > 0, "edgeType"] <- "coexpression"

  # Calculate weights
  combined_interactions$Weight <- rowSums(combined_interactions[, c("experimental", "database", "neighborhood", "cooccurence", "coexpression")])
  combined_interactions$Weight <- combined_interactions$Weight / 1000

  # Create the final edges dataframe from STRINGdb
  combined_edges <- combined_interactions[, c("Gene.1", "Gene.2", "Weight", "edgeType")]

  # Combine STRINGdb and GeneMANIA edges if gm_edges exists
  if(len(db_filepaths) != 0){
    combined_ppi_network <- combined_edges
    for(path in db_filepaths){
      db_edges <- get.DB.edgefile(path)
      combined_ppi_network <- rbind(combined_ppi_network, db_edges)
    }
    assign("ex.ppi_network", combined_ppi_network, envir = .GlobalEnv)} else{ #if db_edges does not exist then do not combine and only use those from STRINGdb
      assign("ex.ppi_network", combined_edges, envir = .GlobalEnv)
    }
}





#' Cluster Filtered Network Example
#' @keywords internal
ex.ClusterFilteredNetwork <- function(accuracy = 0.000001) {
  #Loop through ppi_network and assign every row that matches genenames to an include vector
  include <- c()
  for(a in 1:length(rownames(ex.ppi_network))){
    Gene1 <- ex.ppi_network[a, 1] #Get gene from row a, first col
    Gene2 <- ex.ppi_network[a, 2] #Get gene from row a, second col
    ppi_weight <- as.numeric(ex.ppi_network[a, 3]) #Get the weight of row a

    #If the ppi_weight equals the value that appears in the cccn_matrix for the same genes, add the row number to include
    if(all.equal(ppi_weight, ex.cccn_matrix[Gene1, Gene2], tolerance=accuracy) == TRUE) include[length(include)+1] <- a #all.equal does not return false
  }

  #Assign
  ex.cfn_network <- ex.ppi_network
  if(nrow(ex.cfn_network) == 0) stop("No common edges between PPI edges and cccn_matrix")
  print(ex.cfn_network)
}
