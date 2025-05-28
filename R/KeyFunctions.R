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
  dissimilarity.ptmtable <- 1 - abs(ptmtable.cor)

  # Handle any remaining NA values by setting them to the maximum dissimilarity #
  max_dissimilarity <- max(dissimilarity.ptmtable, na.rm = TRUE)
  dissimilarity.ptmtable[is.na(dissimilarity.ptmtable)] <- max_dissimilarity

  # Make sure the dissimilarity matrix is numeric and suitable for t-SNE #
  dissimilarity.ptmtable <- as.matrix(dissimilarity.ptmtable) #is there a good reason to have this line?

  # Run t-SNE #
  tsne_results <- Rtsne::Rtsne(dissimilarity.ptmtable, dims = 3, perplexity = 15, theta = 0.25, max_iter = 5000, check_duplicates = FALSE, pca = FALSE)
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

  #this is copy and pasted straight from combinedpar so we don't have to run the calculations again -> lmao not anymore
  #no need for its own function I suppose because it's only three lines of code

  #fix spearman thing; so do the exact same thing but no absolute value
  sp.diss.calc <- 1 - ptmtable.cor
  max_diss_sp <- max(sp.diss.calc, na.rm = TRUE)
  sp.diss.calc[is.na(sp.diss.calc)] <- max_dissimilarity
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
}

#' Finds correlations between clusters from multiple distance metrics
#'
#' This function identifies and analyzes clusters using Spearman, Euclidean, and combined
#' t-SNE data, generates cluster size histograms, and saves the plots.
#'
#' @param ptmtable A data frame containing input data for cluster analysis.
#' @param toolong A numeric threshold for cluster separation, defaults to 3.5.
#' @param output_dir The directory where output plots are saved. Defaults to "plots".
#' @export
#'
#' @examples
#' PlotClusters(ptmtable, toolong = 3.5, "output")

PlotClusters <- function(ptmtable, toolong = 3.5, output_dir = "plots") {
    if (!exists("MakeClusterList")) {
        stop("The function 'MakeClusterList' is not defined.")
    }

    # Create output directory if it doesn't exist #
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }

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

}


# Helper function to find intersections of clusters
#'
#' Finds common elements between clusters in three lists.
#'
#' @param list1 A list of ptms.
#' @param list2 A list of ptms.
#' @param list2 A list of ptms.
#' @param keeplength Minimum size of intersections to keep.
#' @return A list of common clusters.
#' @examples
#' FindCommonClusters(list1, list2)
FindCommonClusters <- function(list1, list2, list3 = list(), keeplength = 2){

  #Make sure that the desired column exists in sublists of lists
  if(!TRUE %in% sapply(list1, function(x) "PTM.Name" %in% names(x))) stop("List1 does not have a PTM.Name column")
  if(!TRUE %in% sapply(list2, function(x) "PTM.Name" %in% names(x))) stop("List2 does not have a PTM.Name column")
  if(!TRUE %in% sapply(list3, function(x) "PTM.Name" %in% names(x)) && !!length(list3)) stop("List3 does not have a PTM.Name column") # And list 3 is NOT empty

  #Convert lists into groups of ptms
  list1.ptms <- lapply(list1, function(x){x$"PTM.Name"}) #These are lists of character vectors
  list2.ptms <- lapply(list2, function(y){y$"PTM.Name"})
  list3.ptms <- lapply(list3, function(z){z$"PTM.Name"})

  #Find all the matching intersections of list1 and list2
  returnme <- list()  #Create an empty list to hold those intersections

  for(a in 1:length(list1.ptms)){ #Triple loop to look through elements of the list and compare them
    for(b in 1:length(list2.ptms)){

      #Catch if list3 is empty (for comparing 2 vectors)
      if(!length(list3)){
        temp <- intersect(list1.ptms[[a]], list2.ptms[[b]])
        if(length(temp) > keeplength) returnme[[length(returnme)+1]] <- temp
      }else{
        #if list 3 is NOT empty (for comparing 3 vectors)
        for(c in 1:length(list3.ptms)){
          temp <- Reduce(intersect, list(list1.ptms[[a]], list2.ptms[[b]], list3.ptms[[c]])) #Take the intersection of 3 character vectors (as a vector)
          if(length(temp) > keeplength) returnme[[length(returnme)+1]] <- temp               #And only add it to the list to return if it has enough values
        }
      }
    }
  }
  if(length(returnme) == 0) stop("No common clusters found") #This is for line 370, where the code will return out bounds error anyways if the list is empty!
  return(returnme)
}

#' Replace NAs with 0
#'
#' This function replaces all NA's in a data frame with 0 values.
#'
#' @param df A data frame where NA's are to be replaced with 0s.
#'
#' @return A data frame with NAs replaced by 0.
#' @export
#'
#' @examples
#' NA.to.zero.func(data.frame(a = c(NA, 1), b = c(2, NA)))
NA.to.zero.func <- function(df) {
  df0 <- df
  df0[is.na(df0)] <- 0
  return(df0)
}

#' Generate and Construct All PTMs Network
#'
#' This function generates and constructs the PTMs network from the ptmtable and previously calculated variables.
#'
#' @param ptmtable A data frame containing PTM data.
#' @param keeplength An integer specifying the minimum length of common elements to keep. Default is 2.
#' @param output_dir A string specifying the output directory for saving plots. Default is "plots".
#'
#' @return `eu.sp.sed.ptms`, the list of clusters.
#' @export
#'
#' @examples
#' GenerateAndConstructptmsNetwork(ptmtable)

GenerateAndConstructptmsNetwork <- function(ptmtable, keeplength = 2, output_dir = "plots") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # number of existing entries
  filled   <- function (x)   {length(x) - sum(is.na(x))}

  # Convert lists to data frames #
  #(brackets only take the second and third col and ignore the col that just says number in the cluster) #
  eu.ptms.df  <- plyr::ldply(eu_ptms_list)[, 2:3]
  sp.ptms.df  <- plyr::ldply(sp_ptms_list)[, 2:3]
  sed.ptms.df <- plyr::ldply(sed_ptms_list)[, 2:3]

  # Make group names unique by concatenating a key to the group number (e, s, or sed) #
  eu.ptms.df$group  <- paste(eu.ptms.df$group, "e", sep = "")
  sp.ptms.df$group  <- paste(sp.ptms.df$group, "s", sep = "")
  sed.ptms.df$group <- paste(sed.ptms.df$group, "sed", sep = "")

  # Find common ptms between the three lists using FindCommonClusters()
  eu.sp.sed.ptms <- FindCommonClusters(eu_ptms_list, sp_ptms_list, sed_ptms_list, keeplength)

  #ensure all NAs are zeros
  plot.eu.sp.sed.ptms <- NA.to.zero.func(eu.sp.sed.ptms)

  # Function to generate data frames for heatmaps and evaluations #
  clust.data.from.vec <- function(vec, tbl) {
    if (class(vec) == "list") {
      vec <- unlist(vec)
    }
    at <- tbl[vec, ]
    acol <- names(at[, which(plyr::numcolwise(filled)(at) != 0)])
    if (length(acol) == 1) {
      ats <- data.frame(cbind(rownames(at), as.numeric(at[, acol])))
      names(ats) <- c("PTM.Name", acol)
    } else if (length(acol) >= 2) {
      ats <- cbind(rownames(at), at[, acol])
      names(ats)[1] <- "PTM.Name"
    }
    clust.data <- ats
    return(clust.data)
  } #END

  # Generate data lists for evaluations #
  eu.sp.sed.ptms.data <- list()           #Initilize an empty list
  for (i in 1:length(eu.sp.sed.ptms)) {
    if (length(intersect(eu.sp.sed.ptms[[i]], rownames(ptmtable))) == 0) next  #If the common clusters and rownames of the ptms table have nothing in common break - definitely a better way to write that
    at <- ptmtable[unlist(eu.sp.sed.ptms[[i]]), ] #??? Forgotten paramater?
    if (dim(at)[1] < 2 | dim(at)[2] < 2) next
    eu.sp.sed.ptms.data[[i]] <- clust.data.from.vec(eu.sp.sed.ptms[[i]], tbl = ptmtable)

    # Save the plot
    plot_file <- file.path(output_dir, paste0("plot_", i, ".png"))
    grDevices::png(plot_file, width = 800, height = 600)
    plot(plot.eu.sp.sed.ptms.data[[i]])
    grDevices::dev.off()

    print(paste("Saved plot", i, "to", plot_file))
  }

  # Trim datasets #
  alltrimmedsamples <- apply(ptmtable, 1, filled)
  ptms.t <- ptmtable[which(alltrimmedsamples > 2), ]
  assign('ptmtable', ptms.t)

  return(eu.sp.sed.ptms.data)
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
#' This function orders the adjacency matrix, finds the intersect with the correlation matrix (cccn_matrix), and puts NAs on the diagonal.
#'
#' @param cluster_list A list of clusters to generate adjacency matrices.
#' @param correlation_matrix A correlation matrix to align with the adjacency matrix.
#'
#' @return The CCCN matrix, the intersect between the correlation matrix and the ordered adjacency matrix.
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
  # the adj_matrix appears to be unused in future computations. PERHAPS DELETE?
  return(cccn_matrix)
}
#' Generate Correlation Network
#'
#' This function takes the cccn_matrix, replaces the NA values with zeroes, and turns it into an igraph object.
#'
#' @param cccn_matrix The CCCN matrix, the intersect between the correlation matrix and the ordered adjacency matrix produced from `BindMatrices`.
#'
#' @return An igraph object representing the correlation network.
#' @export
#'
#' @examples
#' CorrelationNetwork(cccn_matrix)
CorrelationNetwork <- function(cccn_matrix) {

  # Make igraph object, replacing NA with 0
  cccn_matrix0 <- NA.to.zero.func(cccn_matrix)
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
#' zero.to.NA.func(data.frame(a = c(NA, 1), b = c(2,NA)))
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
#' @param eu.sp.sed.ptms  A list of all PTMs.
#' @param sed.ptms.peps   A list of SED PTMs peptides.
#' @param AlldataPTMs_cor A correlation matrix for all PTMs.
#'
#' @return A data frame containing PTMs gene correlation edges.
#' @export
#'
#' @examples
#' process_ptms_data(eu.sp.sed.ptms, sed.ptms.peps, AlldataPTMs_cor)
process_ptms_data <- function(eu.sp.sed.ptms, sed.ptms.peps, AlldataPTMs_cor) {
  # Set variables
  eu_sp_sed_ptms <- FindCommonClusters(eu.sp.sed.ptms, sed.ptms.peps, keeplength = 2)

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

#NOT CALLED
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

#' Loads and filters the GeneMania file given a vector of gene names.
#'
#' This helper function loads the GeneMania file
#' and filters out required nodes.
#'
#' @param gmfilename The path to the GeneMania file.
#' @param nodenames A vector containing the names of the relevant genes.
#' @return A data frame with the relevant GeneMania data.
#' @export
#'
#' @examples
#' get.GM.edgefile(gmfilename, nodenames)
get.GM.edgefile <- function(gmfilename, nodenames){

  #reads the file as a table using the first row as a header and tabs as separators (standared for GeneMania interactions)
  gmtable = read.table(gmfilename, header = TRUE, sep = "\t")

  #creates a copy
  gm_edges = gmtable

  #you'll see in about 18 lines
  adjustment = 0

  #iterate through the original table
  for (row in 1:nrow(gmtable)){

    #check to see if both of the genes are in the vector nodenames
    if (gmtable[row, 1] %in% nodenames & gmtable[row, 2] %in% nodenames){

      #nothing happens; I know this is ugly but trust the process

    } else { #if they are NOT in the vector nodenames

      #remove that row (note the adjustment!)
      gm_edges <- gm_edges[-(row + adjustment), ]

      #because we just deleted a row, row 3 in the original is now row 2 for the saved copy
      #so we have to adjust!!
      adjustment = adjustment - 1
    }
  }
  #Removes the column "Network" that just tells what paper this was published in
  gm_edges = gm_edges[ , -5]
  return (gm_edges)
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
  #load("PPIEdges.RData") #Useless?

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

#called :)
# Function to extract gene names from peptide names
pepgene <- function(peps) {
  unique(sapply(peps, function(x) unlist(strsplit(x, " ", fixed=TRUE))[1]))
}

#NOT CALLED
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
genepep.edges.3 <- function(nodelist, pepkey=ld.key) { #NEVER USED
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

