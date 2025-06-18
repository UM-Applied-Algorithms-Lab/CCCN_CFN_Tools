#' Make file for GeneMania input
#'
#' This function outputs a file the user can take to create the GeneMania edgefile.
#'
#' @param cccn_matrix matrix representing the common clusters from the three distance calculations' clusters
#' @param file.path.name path for the output file; defaultsto db_nodes.txt
#'
#' @return A file with all of the gene names which can be copy and pasted into the GeneMania cytoscape app, data frame of the names of the genes
#' @export
#'
#' @examples
#' cccn.cfn.tools:::ex.make_db_input(ex.cccn_matrix)
MakeDBInput <- function(cccn_matrix, file.path.name = "db_nodes.txt") {
  utils::write.table(rownames(cccn_matrix), file = file.path.name, row.names = FALSE, col.names = FALSE, quote = FALSE)
}


#' Process GM Edge File
#'
#' This function process GM edgefile.
#'
#' @param gm.edgefile.path path to GeneMANIA edgefile
#' @param gm.nodetable.path path to GeneMANIA nodetable
#' @param db_nodes.path path to the node file from MakeDBInput
#' @param gm.network.name desired name for the output genemania network; defaulted to gm.network
#'
#' @return GeneMANIA ppi network table
#' @export
#'
#' @examples
#' cccn.cfn.tools::ex.processGMedgefile("path/to/edgefile", "path/to/nodetable", ex.nodenames, ex.gm.network)
ProcessGMEdgefile <- function(gm.edgefile.path, gm.nodetable.path, db_nodes.path, gm.network.name = "gm.network"){
  edgetable <- utils::read.csv(gm.edgefile.path, header = TRUE)        # read the edgefile
  nodetable <- utils::read.csv(gm.nodetable.path, header = TRUE)       # read the nodetable
  nodenames <- utils::read.table(db_nodes.path, header = FALSE)[[1]]   # read the nodenames file

  edgetable <- subset(edgetable, select = c(name, raw.weights))        # only look at the name and raw.weights columns
  nodetable <- subset(nodetable, select = c(name, query.term))         # only look at the names (GM ID) and query.term (real names) columns

  edgetable$Gene.1 <- 'null'   # make new columns!
  edgetable$Gene.2 <- 'null'

  split_names <- strsplit(edgetable$name, "\\|")                                                       # split the names up bc it has three pieces of info there
  edgetable$Gene.1 <- sapply(split_names, function(x)x[1])                                             # take the first thing; first ID name
  edgetable$Gene.1 <- sapply(edgetable$Gene.1, function(x)nodetable$query.term[nodetable$name == x])   # turn the ID into the gene name!
  edgetable$Gene.2 <- sapply(split_names, function(x)x[2])                                             # take the first thing; first ID name
  edgetable$Gene.2 <- sapply(edgetable$Gene.2, function(x)nodetable$query.term[nodetable$name == x])   # turn the ID into the gene name!

  edges <- edgetable[, c("Gene.1", "Gene.2", "raw.weights")]                # sort into a new table with only our information and in the order we want

  colnames(edges) <- c("Gene.1", "Gene.2", "GM.raw.weights")                # rename so when put with everything, clearer where came from

  split_weights <- strsplit(edges$GM.raw.weights, "\\|")                    # split up the weights because they have between 1 and SIX-GODDAMN-TEEN possible weights

  weights.list <- lapply(split_weights, function(x)as.numeric(x))           # make into a number :)
  edges$GM.raw.weights <- sapply(weights.list, function(x)max(x))           # take the maximum? THIS CAN BE CHANGED LATER! I picked max for vibes, but can be mean, sum, etc

  #max_num <- max(edges$GM.raw.weights)                                     # this is code that can be uncommented
  #edges$GM.raw.weights <- sapply(edges$raw.weights, function(x)x/max_num)  # adjusts all of the weights to still have a max of 1; to be used with the sum above ^^^^

  keep <- edges$Gene.1 %in% nodenames & edges$Gene.2 %in% nodenames      # which rows are we keeping
  all.edges <- edges[keep,]                                              # copy 'em over

  assign(gm.network.name, all.edges, envir = .GlobalEnv)                 # assign :)
}
