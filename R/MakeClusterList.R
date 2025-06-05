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
#' MakeClusterList(ptmtable, toolong =  3.5)
MakeClusterList <- function(ptmtable, toolong = 3.5){

  #SPEARMAN CALCULATION

  # Add if statement here to make sure functions are formatted correctly #
  # Ensure ptmtable is a data frame with numeric values #
  ptmtable.sp <- as.data.frame(lapply(ptmtable, as.numeric))

  # Calculate Spearman correlation #
  ptmtable.cor <- stats::cor(t(ptmtable.sp), use = "pairwise.complete.obs", method = "spearman")

  # Replace diagonal with NA #
  diag(ptmtable.cor) <- NA

  # Calculate dissimilarity #
  ptm.correlation.matrix <- 1 - abs(ptmtable.cor)

  # Handle any remaining NA values by setting them to the maximum dissimilarity #
  max_dissimilarity <- max(ptm.correlation.matrix, na.rm = TRUE)
  ptm.correlation.matrix[is.na(ptm.correlation.matrix)] <- max_dissimilarity

  # Make sure the dissimilarity matrix is numeric and suitable for t-SNE #
  colnames(ptm.correlation.matrix) <- rownames(ptmtable) #Repair names
  rownames(ptm.correlation.matrix) <- rownames(ptmtable) #Repair names
  assign("ptm.correlation.matrix", ptm.correlation.matrix, envir = .GlobalEnv) #Correlation Matrix for later use

  # Run t-SNE #
  tsne_results <- Rtsne::Rtsne(ptm.correlation.matrix, dims = 3, perplexity = 15, theta = 0.25, max_iter = 5000, check_duplicates = FALSE, pca = FALSE)
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
  eu.ptms.tsne <- eu.ptms.tsne.list$Y

  # Return the t-SNE results #
  euclidean_result = eu.ptms.tsne

  #COMBINED CALCULATION

  #this is copy and pasted straight from combinedpar so we don't have to run the calculations again
  #no need for its own function I suppose because it's only three lines of code

  #fix spearman thing; so do the exact same thing but no absolute value
  sp.diss.calc <- 1 - ptmtable.cor
  max_diss_sp <- max(sp.diss.calc, na.rm = TRUE)
  sp.diss.calc <- sp.diss.calc * (max_dist / max_diss_sp)
  sp.diss.calc[is.na(sp.diss.calc)] <- 50 * max_dissimilarity
  sp.diss.calc <- as.matrix(sp.diss.calc)

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
    result.span2.df <- data.frame(rownames(ptmtable))
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
  assign("ptm.correlation.matrix", ptmtable.cor, envir = .GlobalEnv)                    #Correlation Matrix for later use
  ########## FIX ME!!!!!! This requires changing the name here back to ptm.cor and then changing every single place it's mentioned in the other files ##########
  ##### Not actually that long but I don't have time rn #####
}
