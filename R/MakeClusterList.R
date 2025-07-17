# Helper function for Rtsne since parameters are set
#
# This function is made to help people new to the Rtsne function understand parameters and organize code. All it does is call Rtsne with set parameters. If you understand R code, feel free to modify the parameters and resource the function.
#
# @param table The matrix that Rtsne needs to be performed on. It's called table because a lot of operations contain the word "matrix" so it was called something else to not be confusing.
# @param iter Number of iterations, defaults 5000
#
# @return A matrix where each row has 3 dimensionial t-SNE coords for a given matrix.
GetRtsne <- function(table, iter=5000){
  if(!is.matrix(table)) table <- as.matrix(table) #Force parameter to be a matrix if it is not already

  # Apply t-SNE to the distance matrix to reduce dimensions to 3 #
  # Parameters: dims = 3 (3D output), perplexity = 15, theta = 0.25 (speed/accuracy trade-off) #
  # max_iter = 5000 (number of iterations)
  # check_duplicates = FALSE (treat rows as unique) #
  # pca = FALSE (no initial PCA) #
  return(Rtsne::Rtsne(table, dims = 3, perplexity = 15, theta = 0.25, max_iter = iter, check_duplicates = FALSE, pca = FALSE))
}

#' Build Clusters Based on PTM Data Under Like Conditions
#'
#' This function takes the post-translational modification table and runs it through three calculations of distance:
#' Euclidean Distance, Spearman Dissimilarity (1 - |Spearman Correlation|), and the average of the two of these.
#' These calculations find the 'distance' between ptms based upon under what conditions they occur.
#' These matricies are then run through t-SNE in order to put them into a 3-dimensional space. Additionally, the intersection of these 3 clusters is also created.
#' A correlation table is also produced based on the Spearman Correlation table.
#'
#' Please note: t-SNE involves an element of randomness; in order to get the same results, set.seed(#) must be called.
#'
#' @param ptmtable A dataset for post-translational modifications. Formatted with numbered rows, and the first column containing PTM names. The rest of the column names should be drugs. Values are numeric values that represent how much the PTM has reacted to the drug.
#' @param keeplength Only keep clusters of ptms whose size is larger than this parameter. (I.e keeplength = 2 then keep ("AARS", "ARMS", "AGRS") but not ("AARS", "ARMS")); default is 2
#' @param correlation.matrix.name Desired name for the correlation matrix to be saved as; defaults to ptm.correlation.matrix
#' @param clusters.list.name Desired name for the lists of clusters to be saved as; defaults to clusters.list
#' @param tsne.coords.name Desired name for the lists of tsne coords to be saved as; defaults to tsne.coords
#' @param common.clusters.name Desired name for the clusters that all 3 methods found in common; defaults to common.clusters
#' @param toolong A numeric threshold for cluster separation, defaults to 3.5.
#' @return The correlation matrix: A data frame showing the correlation between ptms (as the rows and the columns) with NAs placed along the diagonal; and A list of three-dimensional data frames used to represent ptms in space to show relationships between them based on distances. Based on Euclidean Distance, Spearman Dissimilarity, and SED (the average between the two)
#' @export
#'
#' @examples
#' MakeClusterList(ex.ptmtable, "ex.ptm.correlation.matrix", "ex.clusters.list", 3.5)
#' utils::head(ex.ptm.correlation.matrix[, c(1,2,3,4,5)])
#' print(ex.clusters.list[[1]][1])
#' print(ex.clusters.list[[2]][1])
#' print(ex.clusters.list[[3]][1])
MakeClusterList <- function(ptmtable, keeplength = 2, correlation.matrix.name = "ptm.correlation.matrix", clusters.list.name = "clusters.list", tsne.coords.name = "all.tsne.coords", common.clusters.name = "common.clusters", toolong = 3.5){

  #SPEARMAN CALCULATION

  # Add if statement here to make sure functions are formatted correctly #
  # Ensure ptmtable is a data frame with numeric values #
  PTMnames <- sort(rownames(ptmtable))
  ptmtable.sp <- apply(ptmtable, 1:2, as.numeric) #Exclude colnames from matrix and test to make sure they are numbers

  # Calculate Spearman correlation #
  ptm.correlation.matrix <- stats::cor(t(ptmtable.sp), use = "pairwise.complete.obs", method = "spearman")

  # Replace diagonal with NA #
  diag(ptm.correlation.matrix) <- NA

  # Calculate dissimilarity #
  sp.diss.matrix <- 1 - abs(ptm.correlation.matrix)

  # Handle any remaining NA values by setting them to the maximum dissimilarity #
  max.dist.sp <- max(sp.diss.matrix, na.rm = TRUE)

  #Set NA values to 100 * the max distance
  sp.diss.matrix[is.na(sp.diss.matrix)] <- 100 * max.dist.sp

  # Fix names of correlation matrix
  colnames(ptm.correlation.matrix) <- PTMnames
  rownames(ptm.correlation.matrix) <- PTMnames

  # Run t-SNE #
  tsne.results <- GetRtsne(sp.diss.matrix) #Call GetRtsne
  # Return t-SNE results #
  spearman.cluster.coords = tsne.results$Y


  #EUCLIDEAN CALCULATION

  # Add if statement here to make sure functions are formatted correctly #
  # Convert the dataframe to a distance matrix using Euclidean distance #
  ptmtable.dist = as.matrix(stats::dist(ptmtable.sp, method = "euclidean"))

  # Compute the maximum distance in the matrix, excluding NA values #
  max.dist = max(ptmtable.dist, na.rm = TRUE)

  # Replace NA values in the distance matrix with 100 times the maximum distance #
  ptmtable.dist[is.na(ptmtable.dist)] <- 100 * max.dist

  # Normalize the distance matrix by scaling it to a range from 0 to 100. This becomes the distance matrix for euclidian distance which we will run Rtsne on#
  eu.dist.calc <- 100 * ptmtable.dist / max(ptmtable.dist, na.rm = TRUE)
  eu.dist.calc <- as.matrix(eu.dist.calc) #Fix eu.dist.calc RQ

  # Apply t-SNE to the distance matrix to reduce dimensions to 3 #
  # Parameters: dims = 3 (3D output), perplexity = 15, theta = 0.25 (speed/accuracy trade-off) #
  # max_iter = 5000 (number of iterations), check_duplicates = FALSE (treat rows as unique) #
  # pca = FALSE (no initial PCA) #
  eu.ptms.tsne.list <- GetRtsne(eu.dist.calc) #Call GetRtsne

  # Extract the t-SNE results from the output list #
  euclidean.cluster.coords <- eu.ptms.tsne.list$Y

  #COMBINED CALCULATION

  #fix spearman thing; so do the exact same thing but no absolute value
  sp.diss.calc <- 1 - ptm.correlation.matrix            # range goes from (-1 to 1) to (0 to 2)
  max.diss.sp <- max(sp.diss.calc, na.rm = TRUE)        # find the max value (around 2)
  sp.diss.calc[is.na(sp.diss.calc)] <- 50 * max.diss.sp # Set every NA value to 50 * the max dissimilarity
  sp.diss.calc <- as.matrix(sp.diss.calc)               # turn into a matrix

  #find average
  combined.distance <- (sp.diss.calc + eu.dist.calc) / 2
  # Perform t-SNE on the combined distances #
  tsne.result <- GetRtsne(combined.distance) #Call GetRtsne
  sed.cluster.coords <- tsne.result$Y


  #Nested function to analyze result
  clustercreate <- function(result){

    #Compute the minimum spanning tree connecting the points
    tsne.span <- vegan::spantree(stats::dist(result), toolong=toolong)

    #Find clusters that are connected based on toolong (distance?)
    result.disc2 <-  vegan::distconnected(stats::dist(result), toolong = toolong, trace = FALSE)  # test
    # cat ("threshold dissimilarity", toolong, "\n", max(result.disc2), " groups","\n")

    #Create a plot of the clusters using vegan
    vegan::ordiplot(result)
    #lines(tsne.span, result) #???
    vegan::ordihull(result, result.disc2, col="red", lwd=2)

    #Format a data frame
    result.span.df <- data.frame(PTMnames)
    result.span.df$group <- result.disc2 #Add groups found above to the data frame

    #Convert data frame into a list of clusters (check doesn't like group but it's a column name)
    result.span.list <- plyr::dlply(result.span.df, plyr::.(group))  # GROUP LIST  !
    return(result.span.list)

  } #END of nested function

  #Assign tsne coords and distance clusters as a list to global environment
  all.tsne.coords <- list(euclidean.cluster.coords, spearman.cluster.coords, sed.cluster.coords)
  names(all.tsne.coords) <- c("Euclidean", "Spearman", "SED")
  clusters.list <- lapply(all.tsne.coords, clustercreate)
  names(clusters.list) <- c("Euclidean", "Spearman", "SED")


  ### Create common clusters ###
  cluster.as.matrix <- function(cl){ #Helper function to turn a cluster (as a char vector) into a square matrix of 1s
    cl <- cl$PTMnames
    len <- length(cl)
    mat <- matrix(1, len, len)
    rownames(mat) <- cl
    colnames(mat) <- cl
    return(mat)
  }

  #Need Error Catch

  PTMsize <- length(PTMnames)
  clusters.adj.matrix <- matrix(0, PTMsize, PTMsize, dimnames = list(PTMnames, PTMnames)) #Intilize a matirx of 0s

  for(clusters in clusters.list){   #Run for every batch of clusters, run cluster.as.matrix over all of them and turn them into a diagional block matrix

    cl.matrix.list <- lapply(clusters, cluster.as.matrix) #Work on 1 batch of clusters at a time
    temp.matrix <- cl.matrix.list[[1]] #Create the diagonal block matrix
      
    for(i in 2:length(cl.matrix.list)){ #For every cluster in a batch

      addme <- cl.matrix.list[[i]] #This is the matrix we want to add
      row.temp <- nrow(temp.matrix)   #Since the number of rows in temp matrix changes

      temp.matrix <- rbind(temp.matrix, matrix(0, nrow = nrow(addme), ncol = ncol(temp.matrix))) #Make a matrix that is returnmatrix + addme by returnmatrix + addme by filling with 0s
      addme <- rbind(matrix(0, nrow=row.temp, ncol=ncol(addme)), addme)
      temp.matrix <- cbind(temp.matrix, addme)
      rownames(temp.matrix) <- colnames(temp.matrix) #I loathe that I have to do this
    }

    temp.matrix <- temp.matrix[,sort(colnames(temp.matrix))] #Sort the matrices so they line up
    temp.matrix <- temp.matrix[sort(rownames(temp.matrix)),]

    clusters.adj.matrix <- clusters.adj.matrix + temp.matrix #Add all block matrices onto eachother
  }

  clusters.adj.matrix[clusters.adj.matrix < 3] <- NA #Take the combined block matrix and filter out every relationship that is less than 3
  graph <- igraph::graph_from_adjacency_matrix(clusters.adj.matrix) #Create the igraph object
  comp <- igraph::components(graph) #Find the components
  common.clusters <- lapply(seq_along(comp$csize)[comp$csize > keeplength], function(x) igraph::V(graph)$name[comp$membership %in% x]) #Take any greater than size 2 (god bless stack overflow)


  #Assign
  assign(tsne.coords.name, all.tsne.coords, envir = .GlobalEnv) # The list of tsne coords
  assign(clusters.list.name, clusters.list, envir = .GlobalEnv) # list of the t-SNE data for Euclidean, Spearman, and SED
  assign(correlation.matrix.name, ptm.correlation.matrix, envir = .GlobalEnv) # Correlation Matrix for later use
  assign(common.clusters.name, common.clusters, envir = .GlobalEnv) #Common clusters
}
