##WORK IN PROG. JUST COPIED OVER FROM KEYFUNCTIONS##

#' ClusterFilteredNetwork
#'
#' This function finds the intersection between a list of PPI edges and a matrix (which it interprets as a weighted adjacency matrix). Assigns a cfn_network variable to the global environment, which is a list of the common edges between the PPI edges and cccn_matrix.
#' @param accuracy How close a PPI can weight be to a correlation weight to accept the edge. It is recommeneded that genemania users do not increase this parameter as it defaults to the maximum sig figs for a genemania egde weight.
#'
#' @export
#' @examples
#' cccn.cfn.tools:::ex.ClusterFilteredNetwork()
ClusterFilteredNetwork <- function(cccn_matrix, ppi.network, cfn.name = "cfn") {

  #Loop through ppi.network and assign every row that matches genenames to an include vector
  include <- c()
  weights <- c() #CFN
  for(a in 1:length(rownames(ppi.network))){
    #Initilization
    Gene1 <- ppi.network[a, 1] #Get gene from row a, first col
    Gene2 <- ppi.network[a, 2] #Get gene from row a, second col
    if(Gene1 %in% colnames(cccn_matrix) && Gene2 %in% rownames(cccn_matrix)){
      cmw <- cccn_matrix[Gene1, Gene2] #Get the weight
      if(cmw != 0){ #If the weight is nonzero
        include[length(include)+1] <- a
        weights[length(weights)+1] <- cmw
  }}}

  #Assign
  cfn_network <- ppi.network[include, ] #cfn network should only take rows that have been identified by include
  if(nrow(cfn_network) == 0) stop("No common edges between PPI edges and cccn_matrix") #Throw error if no intersection found
  cfn_network$cor_weight <- weights
  assign(cfn.name, cfn_network, envir = .GlobalEnv) #Cluster Filtered Network
}
