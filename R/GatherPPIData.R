#' Make Database Input File
#'
#' This function outputs a file consisting entirely of gene names, each produced on a new line. This data can be copy and pasted into
#' a database input in order to get protein-protein interaction data.
#'
#' @param gene.cccn A matrix showing strength of relationships between proteins using the common clusters between the three distance metrics (Euclidean, Spearman, and Combined (SED))
#' @param file.path.name Path for the output file; defaults to db_nodes.txt
#'
#' @return A file with all of the gene names which can be copy and pasted into the GeneMania cytoscape app, data frame of the names of the genes
#' @export
#'
#' @examples
#' #MakeDBInput(ex.nodenames)
#' cat(ex.nodenames[[1]], sep = '\n')
MakeDBInput <- function(gene.cccn, file.path.name = "db_nodes.txt") {
  utils::write.table(rownames(gene.cccn), file = file.path.name, row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# Pulls nodenames from the gene.cccn
#
# This helper function pulls the gene names from the gene.cccn into a list 'nodenames'
#
# @param gene.cccn A matrix showing strength of relationships between proteins using the common clusters between the three distance metrics (Euclidean, Spearman, and Combined (SED))
# @return data frame of the names of the genes
cccn_to_nodenames <- function(gene.cccn, nodenames.name = 'nodenames'){

  gene.names <- unique(rownames(gene.cccn))

  nodenames <- data.frame(Gene.Names = gene.names, stringsAsFactors = FALSE)

  #return :)
  assign(nodenames.name, nodenames, envir = .GlobalEnv)
  return(nodenames)
}


#' Get STRINGdb PPI data
#'
#' This function finds protein-protein interaction weights by consulting utilizing the STRINGdb database.
#'
#' STRINGdb is required for this function. To download, run:
#' if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#' BiocManager::install("STRINGdb")
#' The full example takes ~10 minutes to load, so it has been commented out and the results are displayed.
#'
#' @param gene.cccn A matrix showing strength of relationships between proteins using common clusters between the three distance metrics (Euclidean, Spearman, and Combined (SED))
#' @param stringdb.name Desired name for the output STRINGdb data frame; defaults to "stringdb.edges"
#' @param nodenames.name Desired name for list of gene names; defaults to nodenames
#'
#' @return Data frame of consisting of the network of interactions from the genes of study pulled from the STRINGdb database and a list of gene names
#' @export
#'
#' @examples
#' # GetSTRINGdb(ex.gene.cccn, 'ex.stringdb.edges', 'ex.nodenames')
#' utils::head(ex.stringdb.edges)
#' utils::head(ex.nodenames)
GetSTRINGdb <- function(gene.cccn, stringdb.name = "stringdb.edges", nodenames.name = "nodenames") {
  nodenames <- cccn_to_nodenames(gene.cccn, nodenames.name)

  if(!requireNamespace("STRINGdb", quietly = TRUE)){
    stop("In order to use this function, please download STRINGdb as described in the vignette, the readme, and the function documentation.")
  }

  # Initialize the STRING database object
  string.db <- STRINGdb::STRINGdb$new(version="12.0", species=9606, score_threshold=0, network_type="full", link_data='full', input_directory="")

  # Retrieve the proteins from the STRING database
  string.proteins <- string.db$get_proteins()

  # Map the genes to STRING IDs
  # please note that nodenames replaces the previous "input_dataset"; nodenames appears to work well :)
  # Gene.Names also replaces experimental
  mapped.genes <- string.db$map(nodenames, "Gene.Names", removeUnmappedRows = TRUE)

  # Retrieve the interactions for the mapped genes
  interactions <- string.db$get_interactions(mapped.genes$STRING_id)

  # Convert protein IDs to gene names
  interactions$Gene.1 <- sapply(interactions$from, function(x) string.proteins[match(x, string.proteins$protein_external_id), "preferred_name"])
  interactions$Gene.2 <- sapply(interactions$to, function(x) string.proteins[match(x, string.proteins$protein_external_id), "preferred_name"])

  # Filter interactions based on evidence types
  str.e <- interactions[interactions$experiments > 0, ]
  str.et <- interactions[interactions$experiments_transferred > 0, ]
  str.d <- interactions[interactions$database > 0, ]
  str.dt <- interactions[interactions$database_transferred > 0, ]

  # Combine filtered interactions
  combined_interactions <- unique(rbind(str.e, str.et, str.d, str.dt))

  # Assign edge types
  combined_interactions$edgeType <- "STRINGdb"
  combined_interactions[combined_interactions$experiments > 0, "edgeType"] <- "experimental"
  combined_interactions[combined_interactions$experiments_transferred > 0, "edgeType"] <- "experimental_transferred"
  combined_interactions[combined_interactions$database > 0, "edgeType"] <- "database"
  combined_interactions[combined_interactions$database_transferred > 0, "edgeType"] <- "database_transferred"

  # Calculate weights
  combined_interactions$Weight <- rowSums(combined_interactions[, c("experiments", "experiments_transferred", "database", "database_transferred")])

  # Create the final edges dataframe from STRINGdb
  combined.edges <- combined_interactions[, c("Gene.1", "Gene.2", "edgeType", "Weight")]
  colnames(combined.edges) <- c("Gene.1", "Gene.2", "Interaction", "STRINGdb.combined_score")

  # assign
  assign(stringdb.name, combined.edges, envir = .GlobalEnv)
}


#' Process GM Edge File
#'
#' This function processes the GM edgefile and translates it back into gene names using the nodetable.
#'
#' @param gm.edgefile.path Path to GeneMANIA edgefile
#' @param gm.nodetable.path Path to GeneMANIA nodetable
#' @param db_nodes.path Path to the node file from MakeDBInput
#' @param gm.network.name Desired name for the output genemania network; defaults to gm.network
#'
#' @return Data frame of consisting of the network of interactions from the genes of study
#' @export
#'
#' @examples
#' ex.edgefile <- system.file("extdata/ex_gm_edgetable.csv", package = "PTMsToPathways")
#' ex.nodefile <- system.file("extdata/ex_gm_nodetable.csv", package = "PTMsToPathways")
#' ex.db.nodes  <- system.file("extdata/ex_db_nodes.txt", package = "PTMsToPathways")
#' ProcessGMEdgefile(ex.edgefile, ex.nodefile, ex.db.nodes, "ex.gm.network")
#' utils::head(ex.gm.network)
ProcessGMEdgefile <- function(gm.edgefile.path, gm.nodetable.path, db_nodes.path, gm.network.name = "gm.network"){
  edgetable <- utils::read.csv(gm.edgefile.path, header = TRUE)        # read the edgefile
  nodetable <- utils::read.csv(gm.nodetable.path, header = TRUE)       # read the nodetable
  nodenames <- utils::read.table(db_nodes.path, header = FALSE)[[1]]   # read the nodenames file

  nodenames <- nodenames[!is.na(nodenames)]   # REMOVE THESE STUPID NAs

  keeper <- edgetable$data.type == "Pathway" | edgetable$data.type == "Physical Interactions"   # which rows have these data types
  edgetable <- edgetable[keeper,]                                                               # copy 'em over

  edgetable <- subset(edgetable, select = c(name, data.type, normalized.max.weight))            # only look at the name and normalized.max.weight columns
  nodetable <- subset(nodetable, select = c(name, query.term))                                  # only look at the names (GM ID) and query.term (real names) columns

  edgetable$Gene.1 <- 'null'   # make new columns!
  edgetable$Gene.2 <- 'null'

  split.names <- strsplit(edgetable$name, "\\|")                                                       # split the names up bc it has three pieces of info there
  edgetable$Gene.1 <- sapply(split.names, function(x)x[1])                                             # take the first thing; first ID name
  edgetable$Gene.1 <- sapply(edgetable$Gene.1, function(x)nodetable$query.term[nodetable$name == x])   # turn the ID into the gene name!
  edgetable$Gene.2 <- sapply(split.names, function(x)x[2])                                             # take the first thing; first ID name
  edgetable$Gene.2 <- sapply(edgetable$Gene.2, function(x)nodetable$query.term[nodetable$name == x])   # turn the ID into the gene name!

  edges <- edgetable[, c("Gene.1", "Gene.2", "data.type", "normalized.max.weight")]      # sort into a new table with only our information and in the order we want

  colnames(edges) <- c("Gene.1", "Gene.2", "interaction", "GM.weights")                  # rename so when put with everything, clearer where came from

  keep <- edges$Gene.1 %in% nodenames & edges$Gene.2 %in% nodenames         # which rows are we keeping
  all.edges <- edges[keep,]                                                 # copy 'em over

  assign(gm.network.name, as.data.frame(all.edges), envir = .GlobalEnv)     # assign :)
}
