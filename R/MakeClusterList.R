#' Populates the global enviroment with cluster groupings based on t-SNE data
#'
#' This function groups t-SNE data points into clusters using a specified threshold
#' and visualizes the clusters.
#'
#' @param ptmtable A dataset for post-translational modifications.
#' @param toolong A numeric threshold for cluster separation, defaults to 3.5.
#' @export
#'
#' @examples
#' ex.MakeClusterList(ex.ptmtable, toolong =  3.5)
MakeClusterList <- function(ptmtable, toolong = 3.5){

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
  assign("eu_ptms_list", clustercreate(euclidean_result, ptmtable), envir = .GlobalEnv) #Matrix containing Euclidean t-SNE coords
  assign("sp_ptms_list", clustercreate(spearman_result, ptmtable), envir = .GlobalEnv)  #Matrix containing Spearman t-SNE coords
  assign("sed_ptms_list", clustercreate(sed_result, ptmtable), envir = .GlobalEnv)      #Matrix containing combined t-SNE coords
  assign("ptm.correlation.matrix", ptm.correlation.matrix, envir = .GlobalEnv)          #Correlation Matrix for later use
}
