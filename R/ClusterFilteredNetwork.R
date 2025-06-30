##WORK IN PROG. JUST COPIED OVER FROM KEYFUNCTIONS##

#' Filter PPI Network by Known Clusters
#'
#' All edges in the PPI network are checked against the cocluster correlation network to ensure existence and a non-zero weight.
#' If either of these conditions are not met, then it will be removed from the list of PPI edges.
#' This new, cluster filtered network is then assigned to the global namespace.
#'
#' @param cccn.matrix A matrix showing strength of relationships between proteins using the common clusters between the three distance metrics (Euclidean, Spearman, and Combined (SED))
#' @param ppi.network A data frame of combined edges from STRINGdb and provided database entries
#' @param cfn.name Desired name of the output cluster filtered network
#'
#' @return cluster filtered network
#'
#' @export
#' @examples
#' ClusterFilteredNetwork(ex.cccn.matrix, ex.ppi.network, cfn.name = "example.cfn")
ClusterFilteredNetwork <- function(cccn.matrix, ppi.network, cfn.name = "cfn") {

  #Loop through ppi.network and assign every row that matches genenames to an include vector
  include <- c()
  weights <- c() #CFN
  for(a in 1:length(rownames(ppi.network))){
    #Initilization
    Gene1 <- ppi.network[a, 1] #Get gene from row a, first col
    Gene2 <- ppi.network[a, 2] #Get gene from row a, second col
    if(Gene1 %in% colnames(cccn.matrix) && Gene2 %in% rownames(cccn.matrix)){
      cmw <- cccn.matrix[Gene1, Gene2] #Get the weight
      if(cmw != 0){ #If the weight is nonzero
        include[length(include)+1] <- a
        weights[length(weights)+1] <- cmw
  }}}

  #Assign
  cfn <- ppi.network[include, ] #cfn should only take rows that have been identified by include
  if(nrow(cfn) == 0) stop("No common edges between PPI edges and cccn.matrix") #Throw error if no intersection found
  cfn$cor_weight <- weights
  assign(cfn.name, cfn, envir = .GlobalEnv) #Cluster Filtered Network
}
