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





#' Process GeneMANIA example
#' @keywords internal
ex.processGMedgefile <- function(epath, npath, ex.nodenames, ex.gm.network){
  gm.ex <- head(ex.gm.network)
  print(gm.ex)
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
