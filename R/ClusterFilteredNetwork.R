##WORK IN PROG. JUST COPIED OVER FROM KEYFUNCTIONS##

#' ClusterFilteredNetwork
#'
#' This function processes correlation edges from a given correlation matrix
#'
#' @return A data frame of correlation edges.
#' @export
#'
#' @examples
#' cccn.cfn.tools:::ClusterFilteredNetwork()
# Function to process correlation edges
ClusterFilteredNetwork <- function() {
  #Error Catch
  if(!exists("cccn_matrix")) stop("cccn_matrix not found. Consider running MakeCorrelationNetwork(keeplength) to make one.")
  if(!exists("ppi_network")) stop("ppi_network not found. Consider looking at step 3 & 4 in https://um-applied-algorithms-lab.github.io/CCCN_CFN_Tools/articles/CoClusterCorrelationNetwork.html to make one.")
  
  #Loop through ppi_network and assign every row that matches genenames to an include vector
  include <- c()
  for(a in 1:length(rownames(ppi_network))){
    Gene1 <- ppi_network[a, 1] #Get gene from row a, first col
    Gene2 <- ppi_network[a, 2] #Get gene from row a, second col
    ppi_weight <- as.numeric(ppi_network[a, 3]) #Get the weight of row a
    
    #If the ppi_weight equals the value that appears in the cccn_matrix for the same genes, add the row number to include
    if(all.equal(ppi_weight, cccn_matrix[Gene1, Gene2], tolerance=.0001) == TRUE) include[length(include)+1] <- a #all.equal does not return false
  }
  
  #Assign
  cfn_network <- ppi_network[include, ]
  assign("cfn_network", cfn_network, envir = .GlobalEnv) #Cluster Filtered Network
}