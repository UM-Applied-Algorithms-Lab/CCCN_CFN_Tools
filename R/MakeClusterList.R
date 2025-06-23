#' Helper function for Rtsne since parameters are set
#' 
#' This function is made to help people new to the Rtsne function understand parameters. 
#' 
#'
#' @param table The matrix that Rtsne needs to be performed on. It's called table because a lot of operations contain the word "matrix" so it was called something else to not be confusing.
#' @param max_iter Number of iterations, defaults 5000
#' 
#' @return A matrix where each row has 3 dimensionial t-SNE coords for a given matrix.  
GetRtsne <- function(table, iter=5000){
  if(!is.matrix(table)) table <- as.matrix(table) #Force parameter to be a matrix if it is not already
  
  # Apply t-SNE to the distance matrix to reduce dimensions to 3 #
  # Parameters: dims = 3 (3D output), perplexity = 15, theta = 0.25 (speed/accuracy trade-off) #
  # max_iter = 5000 (number of iterations)
  # check_duplicates = FALSE (treat rows as unique) #
  # pca = FALSE (no initial PCA) #
  return(Rtsne::Rtsne(table, dims = 3, perplexity = 15, theta = 0.25, max_iter = iter, check_duplicates = FALSE, pca = FALSE)) 
}

#' Populates the global environment with cluster groupings based on t-SNE data
#'
#' This function groups t-SNE data points into clusters using a specified threshold
#' and visualizes the clusters.
#'
#' @param ptmtable A dataset for post-translational modifications. Formatted with numbered rows, and the first column containing PTM names. The rest of the column names should be drugs. Values are numeric values that represent how much the PTM has reacted to the drug. 
#' @param correlation.matrix.name The desired name for the PTM correlation matrix created
#' @param list.name The desired name for the data structure storing the 3 t-SNE matrices
#' @param toolong A numeric threshold for cluster separation, defaults to 3.5.
#' @export
#'
#' @examples
#' MakeClusterList(ex.ptmtable, "example.cor", "example.tsne", 3.5)
MakeClusterList <- function(ptmtable, correlation.matrix.name = "ptm.correlation.matrix", list.name = "clusters.list", toolong = 3.5){

  #SPEARMAN CALCULATION

  # Add if statement here to make sure functions are formatted correctly #
  # Ensure ptmtable is a data frame with numeric values #
  ptmtable.rmnames <- ptmtable[, !colnames(ptmtable) %in% "PTM"] #Remove ANY column that contains strings as such as PTM names, otherwise as.numeric will return NA
  ptmtable.sp <- as.data.frame(lapply(ptmtable.rmnames, as.numeric))

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

  # Run t-SNE #
  tsne_results <- GetRtsne(sp.diss.matrix) #Call GetRtsne
  # Return t-SNE results #
  spearman_cluster_coords = tsne_results$Y


  #EUCLIDEAN CALCULATION

  # Add if statement here to make sure functions are formatted correctly #
  # Convert the dataframe to a distance matrix using Euclidean distance #
  ptmtable.dist = as.matrix(stats::dist(ptmtable.rmnames, method = "euclidean"))

  # Compute the maximum distance in the matrix, excluding NA values #
  max_dist = max(ptmtable.dist, na.rm = TRUE)

  # Replace NA values in the distance matrix with 100 times the maximum distance #
  ptmtable.dist[is.na(ptmtable.dist)] <- 100 * max_dist

  # Normalize the distance matrix by scaling it to a range from 0 to 100. This becomes the distance matrix for euclidian distance which we will run Rtsne on#
  eu.dist.calc <- 100 * ptmtable.dist / max_dist
  eu.dist.calc <- as.matrix(eu.dist.calc) #Fix eu.dist.calc RQ

  # Apply t-SNE to the distance matrix to reduce dimensions to 3 #
  # Parameters: dims = 3 (3D output), perplexity = 15, theta = 0.25 (speed/accuracy trade-off) #
  # max_iter = 5000 (number of iterations), check_duplicates = FALSE (treat rows as unique) #
  # pca = FALSE (no initial PCA) #
  eu.ptms.tsne.list <- GetRtsne(eu.dist.calc) #Call GetRtsne

  # Extract the t-SNE results from the output list #
  euclidean_cluster_coords <- eu.ptms.tsne.list$Y

  #COMBINED CALCULATION

  #fix spearman thing; so do the exact same thing but no absolute value
  sp.diss.calc <- 1 - ptm.correlation.matrix              # range goes from (-1 to 1) to (0 to 2)
  max_diss_sp <- max(sp.diss.calc, na.rm = TRUE)          # find the max value (around 2)
  sp.diss.calc <- sp.diss.calc * (max_dist / max_diss_sp) # SCALING. THIS IS WHERE SCALING OCCURS. All the values are scaled so biggest sp = biggest eu
  sp.diss.calc[is.na(sp.diss.calc)] <- 50 * max_dist_eu   # make the NAs roughly equal to 100
  sp.diss.calc <- as.matrix(sp.diss.calc)                 # turn into a matrix

  #find average
  combined_distance <- (sp.diss.calc + eu.dist.calc) / 2
  # Perform t-SNE on the combined distances #
  tsne_result <- GetRtsne(combined_distance, iter=1000) #Call GetRtsne
  sed_cluster_coords <- tsne_result$Y


  #Nested function to analyze result
  clustercreate <- function(result){

    #Compute the minimum spanning tree connecting the points
    tsne.span <- vegan::spantree(stats::dist(result), toolong=toolong)

    #Find clusters that are connected based on toolong (distance?)
    result.disc2 <-  vegan::distconnected(stats::dist(result), toolong = toolong, trace = TRUE)  # test
    cat ("threshold dissimilarity", toolong, "\n", max(result.disc2), " groups","\n")

    #Create a plot of the clusters using vegan
    vegan::ordiplot(result)
    #lines(tsne.span, result) #???
    vegan::ordihull(result, result.disc2, col="red", lwd=2)

    #Format a data frame
    result.span.df <- data.frame(ptmtable$PTM)
    names(result.span.df) <- "PTM.Name"
    result.span.df$group <- result.disc2 #Add groups found above to the data frame

    #Convert data frame into a list of clusters (check doesn't like group but it's a column name)
    result.span.list <- plyr::dlply(result.span.df, plyr::.(group))  # GROUP LIST  !
    return(result.span.list)

  } #END of nested function

  #Assign distance clusters as a list to global environment
  clusters.list <- list(clustercreate(euclidean_cluster_coords), clustercreate(spearman_cluster_coords), clustercreate(sed_cluster_coords))
  names(clusters.list) <- c("Euclidean", "Spearman", "SED")

  #Assign
  assign(list.name, clusters.list, envir = .GlobalEnv) # list of the t-SNE data for Euclidean, Spearman, and SED
  assign(correlation.matrix.name, ptm.correlation.matrix, envir = .GlobalEnv) # Correlation Matrix for later use
}
