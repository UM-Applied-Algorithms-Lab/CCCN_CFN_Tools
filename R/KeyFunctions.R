load#' Calculates Spearman dissimilarity and t-SNE from a given dataset
#'
#' This function computes the Spearman dissimilarity matrix from the input dataset,
#' processes it for missing values, and performs t-SNE for dimensionality reduction.
#'
#' @param ptmtable A data frame containing numeric data for post-translational modifications.
#' @return A matrix containing t-SNE coordinates (3D).
#' @export
#'
#' @examples
#' SpearmanDissimilarity(ptmtable)
SpearmanDissimilarity <- function(ptmtable) {
    # Add if statement here to make sure functions are formatted correctly #
    # Ensure ptmtable is a data frame with numeric values #
    ptmtable <- as.data.frame(lapply(ptmtable, as.numeric))
    print("Converting Data Types...")

    # Calculate Spearman correlation #
    ptmtable.cor <- stats::cor(t(ptmtable), use = "pairwise.complete.obs", method = "spearman")
    print("Calculating Spearman Correlation...")

    # Replace diagonal with NA #
    diag(ptmtable.cor) <- NA

    # Calculate dissimilarity #
    dissimilarity.ptmtable <- 1 - abs(ptmtable.cor)
    print("Calculating Spearman Dissimilarity...")

    # Handle any remaining NA values by setting them to the maximum dissimilarity #
    max_dissimilarity <- max(dissimilarity.ptmtable, na.rm = TRUE)
    dissimilarity.ptmtable[is.na(dissimilarity.ptmtable)] <- max_dissimilarity
    print("Filtering missing values...")

    # Make sure the dissimilarity matrix is numeric and suitable for t-SNE #
    dissimilarity.ptmtable <- as.matrix(dissimilarity.ptmtable) #is there a good reason to have this line?

    # Run t-SNE #
    tsne_results <- Rtsne::Rtsne(dissimilarity.ptmtable, dims = 3, perplexity = 15, theta = 0.25, max_iter = 5000, check_duplicates = FALSE, pca = FALSE)
    print("Mapping Data Points...")
    # Return t-SNE results #
    return(tsne_results$Y)
}

#' Calculates Euclidean distance and performs t-SNE
#'
#' This function computes the Euclidean distance matrix from the input dataset,
#' normalizes it, and applies t-SNE for dimensionality reduction.
#'
#' @param ptmtable A data frame containing numeric data for post-translational modifications.
#' @return A matrix containing t-SNE coordinates (3D).
#' @export
#'
#' @examples
#' EuclideanDistance(ptmtable)
EuclideanDistance <- function(ptmtable) {
    # Add if statement here to make sure functions are formatted correctly #
    # Convert the dataframe to a distance matrix using Euclidean distance #
    ptmtable.dist = as.matrix(stats::dist(ptmtable, method = "euclidean"))
    print("Converting Data Types...")

    # Compute the maximum distance in the matrix, excluding NA values #
    max_dist = max(ptmtable.dist, na.rm = TRUE)
    print("Finding maximum distance...")

    # Replace NA values in the distance matrix with 100 times the maximum distance #
    ptmtable.dist[is.na(ptmtable.dist)] <- 100 * max_dist
    print("Filtering missing values...")

    # Normalize the distance matrix by scaling it to a range from 0 to 100 #
    ptmtable.dist.1 <- 100 * ptmtable.dist / max_dist
    print("Normalizing distances...")

    # Apply t-SNE to the distance matrix to reduce dimensions to 3 #
    # Parameters: dims = 3 (3D output), perplexity = 15, theta = 0.25 (speed/accuracy trade-off) #
    # max_iter = 5000 (number of iterations), check_duplicates = FALSE (treat rows as unique) #
    # pca = FALSE (no initial PCA) #
    eu.ptms.tsne.list <- Rtsne::Rtsne(as.matrix(ptmtable.dist.1), dims = 3, perplexity = 15, theta = 0.25, max_iter = 5000, check_duplicates = FALSE, pca = FALSE)

    # Extract the t-SNE results from the output list #
    eu.ptms.tsne <- eu.ptms.tsne.list$Y
    print("Mapping Data Points...")

    # Return the t-SNE results #
    return(eu.ptms.tsne)
}

#' Combines Spearman dissimilarity and Euclidean distance in parallel
#'
#' This function uses parallel computing to calculate both Spearman dissimilarity and
#' Euclidean distance, combines them, and performs t-SNE.
#'
#' @param ptmtable A dataset for post-translational modifications.
#' @return A matrix containing t-SNE coordinates (3D).
#' @export
#'
#' @examples
#' CombinedPar(ptmtable)
CombinedPar <- function(ptmtable) {
    # Creates a cluster
    cl <- parallel::makeCluster(2)  # Uses two cores, may increase later #
    # Using makecluster & not parLapply so that this works with Windows machines as well as Unix based ones #
    doParallel::registerDoParallel(cl)

    # Export necessary functions and data to each cluster node #
    parallel::clusterExport(cl, list("SpearmanDissimilarity", "EuclideanDistance", "ptmtable"))
    parallel::clusterEvalQ(cl, {
        library(Rtsne)
        library(parallel)
        library(foreach)
    })

    # Run SpearmanDissimilarity and EuclideanDistance in parallel #
    # Check doesn't like %dopar% and i, also doesn't like foreach::%dopar% -- TODO: figure out.
    results <- foreach::foreach(i = 1:2, .combine = 'list', .packages = c("Rtsne")) %dopar% {
        if (i == 1) {
            return(SpearmanDissimilarity(ptmtable))
        } else {
            return(EuclideanDistance(ptmtable))
        }
    }

    # Extract results #
    spearman_result <- results[[1]]
    euclidean_result <- results[[2]]

    # Continue with the rest of the function #
    combined_distance <- (spearman_result + euclidean_result) / 2

    # Perform t-SNE on the combined distances #
    tsne_result <- Rtsne::Rtsne(as.matrix(combined_distance), dims = 3, perplexity = 15, theta = 0.25, check_duplicates = FALSE, pca = FALSE)
    tsne_coordinates <- tsne_result$Y

    # Stop the cluster #
    parallel::stopCluster(cl)
    return(tsne_coordinates)
}

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

  #find spearman
  spearman_result = SpearmanDissimilarity(ptmtable)
  #find euclidean
  euclidean_result = EuclideanDistance(ptmtable)


  #find average
  #this is copy and pasted straight from combinedpar so we don't have to run the calculations again
  #no need for its own function I suppose because it's only three lines of code
  combined_distance <- (spearman_result + euclidean_result) / 2
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
    names(result.span2.df) <- "Gene.Name"
    result.span2.df$group <- result.disc2 #Add groups found above to the data frame

    #Convert data frame into a list of clusters (check doesn't like group but it's a column name)
    result.span2.list <- plyr::dlply(result.span2.df, plyr::.(group))  # GROUP LIST  !
    return(result.span2.list)

  } #END of nested function

  #Assign different analysises to global enviroment
  assign("eu_ptms_list", clustercreate(euclidean_result, ptmtable), envir = .GlobalEnv) #Matrix containing Euclidean t-SNE coords
  assign("sp_ptms_list", clustercreate(spearman_result, ptmtable), envir = .GlobalEnv)  #Matrix containing Spearman t-SNE coords
  assign("sed_ptms_list", clustercreate(sed_result, ptmtable), envir = .GlobalEnv)      #Matrix containing combined t-SNE coords
}

#' Finds correlations between clusters from multiple distance metrics
#'
#' This function identifies and analyzes clusters using Spearman, Euclidean, and combined
#' t-SNE data, generates cluster size histograms, and saves the plots.
#'
#' @param ptmtable A data frame containing input data for cluster analysis.
#' @param toolong A numeric threshold for cluster separation, defaults to 3.5.
#' @param output_dir The directory where output plots are saved. Defaults to "plots".
#' @return A list containing cluster groupings for each distance metric.
#' @export
#'
#' @examples
#' FindCommonCluster(ptmtable, toolong = 3.5, "output")

FindCommonCluster <- function(ptmtable, toolong = 3.5, output_dir = "plots") {
    if (!exists("MakeClusterList")) {
        stop("The function 'MakeClusterList' is not defined.")
    }

    # Create output directory if it doesn't exist #
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }

    # Make global variables from MakeClusterList #
    MakeClusterList(ptmtable, toolong)

    # Calculate cluster sizes using global variables #
    spsizes_ptms <- sapply(sp_ptms_list, function(x) dim(x)[1])
    sedsizes_ptms <- sapply(sed_ptms_list, function(x) dim(x)[1])
    esizes_ptms <- sapply(eu_ptms_list, function(x) dim(x)[1])

    # Plot and save histograms #
    plot_names <- c("Euclidean_tSNE_Cluster_Sizes.png",
                    "Spearman_tSNE_Cluster_Sizes.png",
                    "Combined_tSNE_Cluster_Sizes.png")

    plot_data <- list(esizes_ptms, spsizes_ptms, sedsizes_ptms)
    plot_colors <- c("yellow", "purple", "brown")
    plot_titles <- c("Euclidean t-SNE Cluster Sizes",
                     "Spearman t-SNE Cluster Sizes",
                     "Combined t-SNE Cluster Sizes")

    for (i in 1:3) {
        grDevices::png(file.path(output_dir, plot_names[i]), width = 800, height = 600)
        graphics::hist(plot_data[[i]], breaks = 100, col = plot_colors[i],
             main = plot_titles[i], xlab = "Cluster Size", ylab = "Frequency")
        grDevices::dev.off()
        print(paste("Saved plot:", plot_names[i]))
    }

    # Return the cluster lists for further use if needed, should also use global variables
    return(list(eu_ptms_list = eu_ptms_list,
                sp_ptms_list = sp_ptms_list,
                sed_ptms_list = sed_ptms_list))
}


# Helper function to find intersections of clusters
#'
#' Finds common elements between clusters in two lists.
#'
#' @param list1 A list of clusters.
#' @param list2 A list of clusters to compare against.
#' @param keeplength Minimum size of intersections to keep.
#' @return A list of common clusters.
#' @examples
#' list.common(cluster_list1, cluster_list2, keeplength = 3)
list.common <- function(list1, list2, keeplength = 3) {
  #Incase of error: 
  #parse <- lapply(list1, function(y) sapply(list2, function(x) intersect(x, y)))
  #prune <- lapply(pare, function(y) return(y[which(sapply(y, function(x) which(length(x) > keeplength)) > 0)]))

  #Parse stores lists of protien clusters. It has lists equal to the SIZE of list1. Each of those lists have clusters equal to the SIZE of list2. The clusters store protiens that are in both list1 and list2.
  parse <- lapply(                 #For a cluster Y, check EVERY cluster in list2 for common protiens and write them to parse. If no protiens in common, list(0)
    list1, function(y){            #Iterate over the clusters in list1 using Y as the iterator
      sapply(list2, function(x){   #Iterate over the clusters in list2 using X as the iterator
        intersect(x, y)            #Returns the common protiens of cluster Y and X. This is done for every cluster in list2 for a cluster Y.  
  })})

  #Similar to the structure of parse, but protien clusters are turned into integers. The integers are equal to the size of the protien cluster in parse.
  dims <- lapply(parse, function(x) sapply(x, length))

  #Identify lists which contain integers greater than keeplength
  keep <- which(sapply(dims, sum) > keeplength)

  #Stores the entire lists which contain valid protien clusters (identified by keep)
  pare <- parse[keep]

  #Since the entire list is stored, filter out all nonvalid clusters in those lists 
  prune <- lapply(
    pare, function(y){                              #For a list Y in pare                       
      return(y[                                     #Find a index of list Y in pare that is a valid cluster and put it in list Y of prune
          which(                                    #Index must be greater than 0 (see the > 0 at the end)
            sapply(y, function(x){                  #For a cluster X in list Y of pare
              which(length(x) > keeplength)}) > 0   #Check if it is has enough protiens to be valid
  )])})
  
  #Create a new list to store all the common clusters? Will return NULL if there is nothing in common.
  newlist <- unlist(prune, recursive = FALSE)
  return(newlist)
}

#' Generate and Construct All PTMs Network
#'
#' This function generates and constructs the PTMs network from given data lists and tables.
#'
#' @param eu_ptms_list A list containing all PTMs data for the Euclidean distance calculation.
#' @param sp_ptms_list A list containing all PTMs data for the Spearman dissimilarity calculation.
#' @param sed_ptms_list A list containing all PTMs data for the average of the Euclidean and Spearman distance calculations.
#' @param ptmtable A data frame containing all PTMs data.
#' @param keeplength An integer specifying the minimum length of common elements to keep. Default is 2.
#' @param output_dir A string specifying the output directory for saving plots. Default is "plots".
#'
#' @return A list containing the updated `ptmtable` and data for `eu.sp.sed.ptms`.
#' @export
#'
#' @examples
#' GenerateAndConstructptmsNetwork(eu_ptms_list, sp_ptms_list, sed_ptms_list, ptmtable)

GenerateAndConstructptmsNetwork <- function(eu_ptms_list, sp_ptms_list, sed_ptms_list,
                                               ptmtable, keeplength = 2, output_dir = "plots") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Mark's Functions #
  without  <- function(x, y)  x[!x %in% y]             # x without y
  nmissing <- function(x)     sum(is.na(x))            # number of missing entries
  filled   <- function (x)   {length(x) - nmissing(x)} # number of existing entries
  fractNA  <- function(df) {
      result <- nmissing(df)/(dim(df)[1]*dim(df)[2])  #fraction of total entries that are missing
      return(result)
  }
  mean.na   <- function(x)     mean(x, na.rm=TRUE) #why do these exist? 
  max.na    <- function(x)     max(x, na.rm=TRUE)
  min.na    <- function(x)     min(x, na.rm=TRUE)
  sd.na     <- function(x)     stats::sd(x, na.rm=TRUE)
  outersect <- function(x,y)  {sort(c(setdiff(x,y), setdiff(y,x)))}

  # Convert lists to data frames #
  #(brackets only take the second and third col and ignore the col that just says number in the cluster) #
  eu.ptms.df  <- plyr::ldply(eu_ptms_list)[, 2:3]
  sp.ptms.df  <- plyr::ldply(sp_ptms_list)[, 2:3]
  sed.ptms.df <- plyr::ldply(sed_ptms_list)[, 2:3]

  # Make group names unique by concatenating a key to the group number (e, s, or sed) #
  eu.ptms.df$group  <- paste(eu.ptms.df$group, "e", sep = "")
  sp.ptms.df$group  <- paste(sp.ptms.df$group, "s", sep = "")
  sed.ptms.df$group <- paste(sed.ptms.df$group, "sed", sep = "")

  # Group everything together, data frames pasted together with rows of E on top of rows of Sed on top of rows of S #
  ptmsgroups.df <- rbind(eu.ptms.df, sed.ptms.df, sp.ptms.df)

  # Functions to extract gene names and PTMs #
  extract.genes.from.clist <- function(clusterlist.element) {
    element <- clusterlist.element[1]
    genes   <- unique(sapply(as.character(element$Gene.Name), function(x) unlist(strsplit(x, " ", fixed = TRUE))[1]))
    return(genes)
  }

  extract.peps.from.clist <- function(clusterlist.element) {
    element <- clusterlist.element[1]
    return(as.character(element$Gene.Name))
  }

  eu.ptms.genes  <- extract.genes.from.clist(eu.ptms.df)
  sp.ptms.genes  <- extract.genes.from.clist(sp.ptms.df)
  sed.ptms.genes <- extract.genes.from.clist(sed.ptms.df)

  eu.ptms.peps  <- extract.peps.from.clist(eu.ptms.df)
  sp.ptms.peps  <- extract.peps.from.clist(eu.ptms.df)
  sed.ptms.peps <- extract.peps.from.clist(eu.ptms.df)


  ################# FAILING HERE #################


  eu.sp.ptms <- list.common(eu.ptms.peps, sp.ptms.peps, keeplength)
  eu.sp.ptms.sizes <- sapply(eu.sp.ptms, length)
  eu.sp.sed.ptms <- list.common(eu.sp.ptms, sed.ptms.peps, keeplength) #nonnull when keeplength < 2
  eu.sp.sed.ptms.sizes <- sapply(eu.sp.sed.ptms, length)


  ################# FAILING HERE #################


  # Function to generate data frames for heatmaps and evaluations #
  clust.data.from.vec <- function(vec, tbl) {
    if (class(vec) == "list") {
      vec <- unlist(vec)
    }
    at <- tbl[vec, ]
    acol <- names(at[, which(plyr::numcolwise(filled)(at) != 0)])
    if (length(acol) == 1) {
      ats <- data.frame(cbind(rownames(at), as.numeric(at[, acol])))
      names(ats) <- c("Gene.Name", acol)
    } else if (length(acol) >= 2) {
      ats <- cbind(rownames(at), at[, acol])
      names(ats)[1] <- "Gene.Name"
    }
    clust.data <- ats
    return(clust.data)
  }

  # Generate data lists for evaluations #
  eu.sp.sed.ptms.data <- list()
  for (i in 1:length(eu.sp.sed.ptms)) {
    if (length(intersect(eu.sp.sed.ptms[[i]], rownames(ptmtable))) == 0) next
    at <- ptmtable[unlist(eu.sp.sed.ptms[[i]]), ]
    if (dim(at)[1] < 2 | dim(at)[2] < 2) next
    eu.sp.sed.ptms.data[[i]] <- clust.data.from.vec(eu.sp.sed.ptms[[i]], tbl = ptmtable)

    # Save the plot
    plot_file <- file.path(output_dir, paste0("plot_", i, ".png"))
    grDevices::png(plot_file, width = 800, height = 600)
    plot(eu.sp.sed.ptms.data[[i]])
    grDevices::dev.off()

    print(paste("Saved plot", i, "to", plot_file))
  }

  # Trim datasets #
  alltrimmedsamples <- apply(ptmtable, 1, filled)
  ptms.t <- ptmtable[which(alltrimmedsamples > 2), ]
  ptmtable <- ptms.t

  # Repair bad clusters #
  bad.clusterlist <- list()
  badptms <- unique(outersect(rownames(ptmtable), rownames(ptmtable)))

  return(list(ptmtable = ptmtable, eu.sp.sed.ptms.data = eu.sp.sed.ptms.data))
}

#' Create Adjacency Matrix
#'
#' This function creates an adjacency matrix for a given list element.
#'
#' @param list.element A list of elements to construct the adjacency matrix.
#'
#' @return A square matrix where rows and columns correspond to the input list elements.
#' @export
#'
#' @examples
#' MakeAdjMatrix(c("A", "B", "C"))
MakeAdjMatrix <- function(list.element) {
  list.el.mat <- matrix(1, nrow = length(list.element), ncol = length(list.element))
  rownames(list.el.mat) <- list.element
  colnames(list.el.mat) <- list.element
  return(list.el.mat)
}
#' Bind Matrices
#'
#' This function binds matrices, aligns them, and prepares adjacency and CCCN matrices.
#'
#' @param cluster_list A list of clusters to generate adjacency matrices.
#' @param correlation_matrix A correlation matrix to align with the adjacency matrix.
#'
#' @return A list containing the combined adjacency matrix and CCCN matrix.
#' @export
#'
#' @examples
#' BindMatrices(cluster_list, correlation_matrix)
BindMatrices <- function(cluster_list, correlation_matrix) {
  # Generate the combined adjacency matrix
  adj_matrix <- plyr::rbind.fill.matrix(plyr::llply(cluster_list, MakeAdjMatrix))
  rownames(adj_matrix) <- colnames(adj_matrix)

  # Order the adjacency matrix by row and column names
  adj_matrix_ordered <- adj_matrix[order(rownames(adj_matrix)), order(colnames(adj_matrix))]

  # Align the correlation matrix with the ordered adjacency matrix
  matched_rows <- intersect(rownames(adj_matrix_ordered), rownames(correlation_matrix))
  matched_cols <- intersect(colnames(adj_matrix_ordered), colnames(correlation_matrix))
  cccn_matrix <- correlation_matrix[matched_rows, matched_cols]

  # Replace NA values in the correlation matrix
  na_indices <- which(is.na(adj_matrix_ordered), arr.ind = TRUE)
  cccn_matrix <- replace(cccn_matrix, na_indices, NA)

  # Remove self-loops by setting diagonal to NA
  if (any(!is.na(diag(cccn_matrix)))) {
    diag(cccn_matrix) <- NA
  }

  # Return the adjacency and CCCN matrices as a list
  return(list(adj_matrix = adj_matrix_ordered, cccn_matrix = cccn_matrix))
}
#' Generate Correlation Network
#'
#' This function creates a correlation network graph from a given set of matrices.
#'
#' @param bind_result A list containing the adjacency and CCCN matrices from `BindMatrices`.
#'
#' @return An igraph object representing the correlation network.
#' @export
#'
#' @examples
#' CorrelationNetwork(bind_result)
CorrelationNetwork <- function(bind_result) {
  adj_matrix <- bind_result$adj_matrix
  cccn_matrix <- bind_result$cccn_matrix

  # Make igraph object, replacing NA with 0
  cccn_matrix0 <- cccn_matrix
  cccn_matrix0[is.na(cccn_matrix0)] <- 0
  graph <- igraph::graph_from_adjacency_matrix(as.matrix(cccn_matrix0), mode = "lower", diag = FALSE, weighted = "Weight")

  # Return the graph object
  return(graph)
}

#' Replace Zeros with NA
#'
#' This function replaces all zeros in a data frame with NA values.
#'
#' @param df A data frame where zeros are to be replaced with NA.
#'
#' @return A data frame with zeros replaced by NA.
#' @export
#'
#' @examples
#' zero.to.NA.func(data.frame(a = c(0, 1), b = c(2, 0)))
zero.to.NA.func <- function(df) {
  cf <- df
  zer0 <- which(cf==0, arr.ind = TRUE)
  cfNA <- as.matrix(cf)
  cfNA[zer0] <- NA
  cfNA <- data.frame(cfNA)
  return(cfNA)
}

#' Process PTMs Data
#'
#' This function processes PTMs data, creates correlation networks, and constructs adjacency matrices.
#'
#' @param eu.sp.sed.ptms A list of all PTMs.
#' @param sed.ptms.peps A list of SED PTMs peptides.
#' @param AlldataPTMs_cor A correlation matrix for all PTMs.
#'
#' @return A data frame containing PTMs gene correlation edges.
#' @export
#'
#' @examples
#' process_ptms_data(eu.sp.sed.ptms, sed.ptms.peps, AlldataPTMs_cor)
process_ptms_data <- function(eu.sp.sed.ptms, sed.ptms.peps, AlldataPTMs_cor) {
  # Set variables
  eu_sp_sed_ptms <- list.common(eu.sp.sed.ptms, sed.ptms.peps, keeplength = 2)

  # Create adjacency matrices
  ptms_adj <- plyr::rbind.fill.matrix(plyr::llply(eu_sp_sed_ptms, MakeAdjMatrix))
  rownames(ptms_adj) <- colnames(ptms_adj)

  # Order and align matrices
  ptms_adj_o <- ptms_adj[order(rownames(ptms_adj)), order(colnames(ptms_adj))]

  ptms_cccn_1 <- AlldataPTMs_cor[rownames(AlldataPTMs_cor) %in% rownames(ptms_adj_o), colnames(AlldataPTMs_cor) %in% colnames(ptms_adj_o)]

  # Check matrices
  if(length(setdiff(rownames(ptms_adj), rownames(ptms_cccn_1))) != 0) stop("Mismatch in rownames")
  if(length(intersect(rownames(ptms_adj), rownames(AlldataPTMs_cor))) != nrow(ptms_adj)) stop("Mismatch in intersect rownames")

  # Add correlation as edge values in adjacency matrix
  ptms_cccn <- AlldataPTMs_cor[intersect(rownames(ptms_adj_o), rownames(AlldataPTMs_cor)), intersect(colnames(ptms_adj_o), colnames(AlldataPTMs_cor))]

  # Replace NA values
  ptms_NA <- which(is.na(ptms_adj_o), arr.ind = TRUE)
  ptms_cccn <- replace(ptms_cccn, ptms_NA, NA)
  if (any(!is.na(diag(ptms_cccn)))) diag(ptms_cccn) <- NA

  # Make igraph objects
  ptms_cccn0 <- ptms_cccn
  ptms_cccn0[is.na(ptms_cccn0)] <- 0
  ptms_cccn_g <- igraph::graph_from_adjacency_matrix(as.matrix(ptms_cccn0), mode = "lower", diag = FALSE, weighted = "Weight")

  # Gene CCCN construction
  ptms_gene_cccn <- data.frame(ptms_cccn, row.names = rownames(ptms_cccn), check.rows = TRUE, check.names = FALSE, fix.empty.names = FALSE)
  ptms_gene_cccn$Gene_Name <- sapply(rownames(ptms_gene_cccn), function(x) unlist(strsplit(x, " ", fixed = TRUE))[1])

  ptms_gene_cccn[lower.tri(ptms_gene_cccn)] <- NA

  #check doesn't like Gene_Name but it's a column name
  ptms_gene_cccn2 <- plyr::ddply(ptms_gene_cccn, plyr::.(Gene_Name), plyr::numcolwise(function(x) sum(x, na.rm = TRUE)), .progress = "tk")

  rownames(ptms_gene_cccn2) <- ptms_gene_cccn2$Gene_Name
  ptms_gene_cccn2 <- ptms_gene_cccn2[, 2:ncol(ptms_gene_cccn2)]
  ptms_gene_cccn2 <- data.frame(t(ptms_gene_cccn2))
  ptms_gene_cccn2$Gene <- sapply(rownames(ptms_gene_cccn2), function(x) unlist(strsplit(x, " ", fixed = TRUE))[1])

  #check doesn't like Gene but it's a column name
  ptms_gene_cccn3 <- plyr::ddply(ptms_gene_cccn2, plyr::.(Gene), plyr::numcolwise(function(x) sum(x, na.rm = TRUE)), .progress = "tk")

  names(ptms_gene_cccn3)[2:ncol(ptms_gene_cccn3)] <- ptms_gene_cccn3$Gene
  rownames(ptms_gene_cccn3) <- ptms_gene_cccn3$Gene

  ptms_gene_cccn0 <- ptms_gene_cccn3[, 2:ncol(ptms_gene_cccn3)]
  ptms_gene_cccn_na <- zero.to.NA.func(ptms_gene_cccn0)

  ptms_gene_cccn_g <- igraph::graph.adjacency(as.matrix(ptms_gene_cccn0), mode = "lower", diag = FALSE, weighted = "Weight")

  ptms_gene_cccn_edges <- data.frame(igraph::as_edgelist(ptms_gene_cccn_g))
  names(ptms_gene_cccn_edges) <- c("Gene.1", "Gene.2")
  ptms_gene_cccn_edges$Weight <- igraph::edge_attr(ptms_gene_cccn_g)[[1]]
  ptms_gene_cccn_edges$interaction <- "correlation"
  ptms_gene_cccn_edges$interaction[ptms_gene_cccn_edges$Weight <= -0.5] <- "negative correlation"
  ptms_gene_cccn_edges$interaction[ptms_gene_cccn_edges$Weight >= 0.5] <- "positive correlation"

  return(ptms_gene_cccn_edges)
}

#' Extract Gene Names from Peptide Vector
#'
#' This function extracts gene names from a given peptide vector.
#'
#' @param pepvec A vector of peptides.
#' @param pepsep A string specifying the separator for peptides. Default is "; ".
#'
#' @return A vector of unique gene names.
#' @export
#'
#' @examples
#' get.gene.names.from.peps(c("gene1 peptide1", "gene2 peptide2"))
get.gene.names.from.peps <- function(pepvec, pepsep="; ") {
  genevec=NULL
  for(i in 1:length(pepvec)) {
    x <- unlist(strsplit(as.character(pepvec[i]), pepsep, fixed=TRUE))
    genes <- unique(sapply(as.character(x),  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1]))
    genevec <- c(genevec, genes)
  }
  return(genevec)
}

#' Find PPI Edges
#'
#' This function finds protein-protein interaction edges by combining input datasets with STRING and GeneMANIA databases.
#'
#' @param input_dataset The input dataset containing experimental data.
#' @param gmfilename The filename of the GeneMANIA data.
#' @param nodenames A vector of node names.
#'
#' @return A data frame of combined edges from STRINGdb and GeneMANIA.
#' @export
#'
#' @examples
#' find_ppi_edges("input_data.txt", "gmfilename.txt", nodenames)
find_ppi_edges <- function(input_dataset, gmfilename, nodenames) {
  # Load PPI edges from other databases
  load("PPIEdges.RData")

  # Initialize the STRING database object
  string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=0, link_data="detailed", input_directory="")

  # Retrieve the proteins from the STRING database
  string_proteins <- string_db$get_proteins()
  print(dim(string_proteins))

  # Read the dataset that you want to combine with the STRING database
  filter_db <- utils::read.table(input_dataset, header = TRUE, sep = "\t")
  print(colnames(filter_db))

  if (!"experimental" %in% colnames(filter_db)) {
    stop("Column 'experimental' not found in input dataset.")
  }

  # Map the genes to STRING IDs
  mapped_genes <- string_db$map(filter_db, "experimental", removeUnmappedRows = TRUE)
  print(utils::head(mapped_genes))

  # Retrieve the interactions for the mapped genes
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)

  # Convert protein IDs to gene names
  interactions$Gene.1 <- sapply(interactions$from, function(x) string_proteins[match(x, string_proteins$protein_external_id), "preferred_name"])
  interactions$Gene.2 <- sapply(interactions$to, function(x) string_proteins[match(x, string_proteins$protein_external_id), "preferred_name"])

  # Filter interactions based on evidence types
  str.e <- interactions[interactions$experiments > 0, ]
  str.et <- interactions[interactions$experiments_transferred > 0, ]
  str.d <- interactions[interactions$database > 0, ]
  str.dt <- interactions[interactions$database_transferred > 0, ]

  # Combine filtered interactions
  combined_interactions <- unique(rbind(str.e, str.et, str.d, str.dt))

  # Assign edge types
  combined_interactions$edgeType <- "STRINGdb"
  combined_interactions[combined_interactions$database > 0, "edgeType"] <- "database"
  combined_interactions[combined_interactions$database_transferred > 0, "edgeType"] <- "database"
  combined_interactions[combined_interactions$experiments > 0, "edgeType"] <- "experiments"
  combined_interactions[combined_interactions$experiments_transferred > 0, "edgeType"] <- "experiments"

  # Calculate weights
  combined_interactions$Weight <- rowSums(combined_interactions[, c("experiments", "experiments_transferred", "database", "database_transferred")])
  combined_interactions$Weight <- combined_interactions$Weight / 1000

  # Create the final edges dataframe from STRINGdb
  combined_edges <- combined_interactions[, c("Gene.1", "Gene.2", "Weight", "edgeType")]

  # Get GeneMANIA edges
  gm_edges <- get.GM.edgefile(gmfilename, nodenames)

  # Combine STRINGdb and GeneMANIA edges
  final_edges <- rbind(combined_edges, gm_edges)

  return(final_edges)
}

# Function to extract gene names from peptide names
pepgene <- function(peps) {
  unique(sapply(peps, function(x) unlist(strsplit(x, " ", fixed=TRUE))[1]))
}

#' Extract Gene Names from Peptide Edge File
#'
#' This function extracts unique gene names from a peptide edge file.
#'
#' @param peptide.edgefile A data frame containing peptide edge information.
#'
#' @return A vector of unique gene names.
#' @export
#'
#' @examples
#' extract.gene.names(peptide.edgefile)
# Function to extract gene names from peptide edge file
extract.gene.names <- function(peptide.edgefile) {
  peps <- c(peptide.edgefile[,1], peptide.edgefile[,2])
  genes <- unique(sapply(peps, function(x) unlist(strsplit(x, " ", fixed=TRUE))[1]))
  return(genes)
}

#' Create Gene-Peptide Edges
#'
#' This function creates peptide edges for a given node list.
#'
#' @param nodelist A vector of node names.
#' @param pepkey A data frame containing peptide keys.
#'
#' @return A data frame of peptide edges with weights and edge types.
#' @export
#'
#' @examples
#' genepep.edges.3(nodelist, pepkey)
genepep.edges.3 <- function(nodelist, pepkey=ld.key) {
  nodelist <- unique(nodelist)
  gpedges <- pepkey[pepkey$Gene.Name %in% nodelist, 1:2]
  names(gpedges)[1:2] <- c("Gene.1", "Gene.2")
  gpedges$edgeType <- "peptide"
  gpedges$Weight <- 1
  gpedges$Alt.Weight <- 100
  gpedges$Directed <- FALSE
  return(unique(gpedges))
}

#' Process Correlation Edges
#'
#' This function processes correlation edges from a given correlation matrix.
#'
#' @param cor_matrix A correlation matrix.
#' @param mode A string specifying the graph mode. Default is "lower".
#'
#' @return A data frame of correlation edges.
#' @export
#'
#' @examples
#' process_correlation_edges(cor_matrix)
# Function to process correlation edges
process_correlation_edges <- function(cor_matrix, mode="lower") {
  g <- igraph::graph_from_adjacency_matrix(as.matrix(cor_matrix), mode=mode, diag=FALSE, weighted="Weight")
  edges <- data.frame(igraph::as_edgelist(g))
  edges$Weight <- igraph::edge_attr(g)[[1]]
  edges$edgeType <- "correlation"
  edges$edgeType[edges$Weight <= -0.5] <- "negative correlation"
  edges$edgeType[edges$Weight >= 0.5] <- "positive correlation"
  edges <- edges[!is.na(edges$Weight),]
  names(edges)[1:2] <- c("Peptide.1", "Peptide.2")
  edges$Gene.1 <- sapply(edges$Peptide.1, pepgene)
  edges$Gene.2 <- sapply(edges$Peptide.2, pepgene)
  return(edges)
}

# Function to filter dual modifications
filter_dual_modifications <- function(edges, mod1, mod2) {
  dual_mod <- edges[intersect(grep(mod1, edges$Peptide.1), grep(mod2, edges$Peptide.2)), ]
  return(dual_mod)
}

# Function to analyze negative correlations
analyze_negative_correlations <- function(edges) {
  neg <- edges[edges$Weight < 0, ]
  vneg <- neg[abs(neg$Weight) >= 0.5, ]
  vvneg <- neg[abs(neg$Weight) > 0.543, ]

  neg_genes <- unique(neg$Gene.1)
  vneg_genes <- unique(vneg$Gene.1)
  vvneg_genes <- unique(vvneg$Gene.1)

  return(list(neg=neg, vneg=vneg, vvneg=vvneg,
              neg_genes=neg_genes, vneg_genes=vneg_genes, vvneg_genes=vvneg_genes))
}

