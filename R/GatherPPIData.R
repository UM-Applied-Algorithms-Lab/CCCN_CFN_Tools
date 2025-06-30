#' Make Database Input File
#'
#' This function outputs a file consisting entirely of gene names, each produced on a new line. This data can be copy and pasted into
#' a database input in order to get protein-protein interaction data.
#'
#' @param cccn.matrix Matrix representing the common clusters from the three distance calculations' clusters
#' @param file.path.name Path for the output file; defaults to db_nodes.txt
#'
#' @return A file with all of the gene names which can be copy and pasted into the GeneMania cytoscape app, data frame of the names of the genes
#' @export
#'
#' @examples
#' cccn.cfn.tools:::ex.MakeDBInput(ex.nodenames)
MakeDBInput <- function(cccn.matrix, file.path.name = "db_nodes.txt") {
  utils::write.table(rownames(cccn.matrix), file = file.path.name, row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# Pulls nodenames from the cccn.matrix
#
# This helper function pulls the gene names from the cccn.matrix into a list 'nodenames'
#
# @param cccn.matrix dataframe of dataframes that represent the common clusters from the three distance calculations' clusters
# @return data frame of the names of the genes
cccn_to_nodenames <- function(cccn.matrix, nodenames.name = 'nodenames'){

  gene.names <- unique(rownames(cccn.matrix))

  nodenames <- data.frame(Gene.Names = gene.names, stringsAsFactors = FALSE)

  #return :)
  assign(nodenames.name, nodenames, envir = .GlobalEnv)
}


#' Get STRINGdb PPI data
#'
#' This function finds protein-protein interaction weights by consulting utilizing the STRINGdb database.
#'
#' @param cccn.matrix A matrix showing strength of relationships between proteins using common clusters between the three distance metrics (Euclidean, Spearman, and Combined (SED))
#' @param STRINGdb.name Desired name for the output STRINGdb data frame; defaults to "string.edges"
#' @param nodes.name Desired name for list of gene names; defaults to nodenames
#'
#' @return Data frame of consisting of the network of interactions from the genes of study and a list of gene names
#' @export
#'
#' @examples
#' cccn.cfn.tools:::ex.FindPPIEdges(ex.cccn.matrix, "ex.string.edges", "ex.nodenames")
GetSTRINGdb <- function(cccn.matrix, STRINGdb.name = "string.edges", nodenames.name = "nodenames") {
  cccn_to_nodenames(cccn.matrix, nodenames.name)

  if (!exists("STRINGdb")){                          # check if stringdb is libraried

    if(system.file(package="STRINGdb") == ""){       # check if stringdb is installed at all

      BiocManager::install("STRINGdb")               # install!
    }
    library(STRINGdb)                                # library it
  }

  # Initialize the STRING database object
  string.db <- STRINGdb$new(version="12.0", species=9606, score_threshold=0, input_directory="")

  # Retrieve the proteins from the STRING database
  string.proteins <- string.db$get_proteins()

  if (!"Gene.Names" %in% colnames(nodenames)) {
    stop("Column 'Gene.Names' not found in nodenames.")
  }

  # Map the genes to STRING IDs
  # please note that nodenames replaces the previous "input_dataset"; nodenames appears to work well :)
  # Gene.Names also replaces experimental
  mapped.genes <- string.db$map(nodenames, "Gene.Names", removeUnmappedRows = TRUE)
  print(utils::head(mapped.genes))

  # Retrieve the interactions for the mapped genes
  interactions <- string.db$get_interactions(mapped.genes$STRING_id)

  # Convert protein IDs to gene names
  interactions$Gene.1 <- sapply(interactions$from, function(x) string.proteins[match(x, string.proteins$protein_external_id), "preferred_name"])
  interactions$Gene.2 <- sapply(interactions$to, function(x) string.proteins[match(x, string.proteins$protein_external_id), "preferred_name"])

  # Create the final edges dataframe from STRINGdb
  combined.edges <- interactions[, c("Gene.1", "Gene.2", "combined_score")]
  colnames(combined.edges) <- c("Gene.1", "Gene.2", "STRINGdb.combined_score")

  # assign
  assign(STRINGdb.name, combined.edges, envir = .GlobalEnv)
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
#' cccn.cfn.tools:::ex.ProcessGMEdgefile("edgefile.csv", "nodetable.csv", ex.nodenames, ex.gm.network)
ProcessGMEdgefile <- function(gm.edgefile.path, gm.nodetable.path, db_nodes.path, gm.network.name = "gm.network"){
  edgetable <- utils::read.csv(gm.edgefile.path, header = TRUE)        # read the edgefile
  nodetable <- utils::read.csv(gm.nodetable.path, header = TRUE)       # read the nodetable
  nodenames <- utils::read.table(db_nodes.path, header = FALSE)[[1]]   # read the nodenames file

  nodenames <- nodenames[!is.na(nodenames)]   # REMOVE THESE STUPID NAs

  keeper <- edgetable$data.type == "Pathway" | edgetable$data.type == "Physical Interactions"   # which rows have these data types
  edgetable <- edgetable[keeper,]                                                               # copy 'em over

  edgetable <- subset(edgetable, select = c(name, normalized.max.weight))        # only look at the name and normalized.max.weight columns
  nodetable <- subset(nodetable, select = c(name, query.term))         # only look at the names (GM ID) and query.term (real names) columns

  edgetable$Gene.1 <- 'null'   # make new columns!
  edgetable$Gene.2 <- 'null'

  split.names <- strsplit(edgetable$name, "\\|")                                                       # split the names up bc it has three pieces of info there
  edgetable$Gene.1 <- sapply(split.names, function(x)x[1])                                             # take the first thing; first ID name
  edgetable$Gene.1 <- sapply(edgetable$Gene.1, function(x)nodetable$query.term[nodetable$name == x])   # turn the ID into the gene name!
  edgetable$Gene.2 <- sapply(split.names, function(x)x[2])                                             # take the first thing; first ID name
  edgetable$Gene.2 <- sapply(edgetable$Gene.2, function(x)nodetable$query.term[nodetable$name == x])   # turn the ID into the gene name!

  edges <- edgetable[, c("Gene.1", "Gene.2", "normalized.max.weight")]      # sort into a new table with only our information and in the order we want

  colnames(edges) <- c("Gene.1", "Gene.2", "GM.weights")                    # rename so when put with everything, clearer where came from

  keep <- edges$Gene.1 %in% nodenames & edges$Gene.2 %in% nodenames         # which rows are we keeping
  all.edges <- edges[keep,]                                                 # copy 'em over

  assign(gm.network.name, as.data.frame(all.edges), envir = .GlobalEnv)                 # assign :)
}
