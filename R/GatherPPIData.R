#' Make file for GeneMania input
#'
#' This function outputs a file the user can take to create the GeneMania edgefile.
#'
#' @param cccn_matrix matrix representing the common clusters from the three distance calculations' clusters
#' @param file.path.name path for the output file; defaults to db_nodes.txt
#'
#' @return A file with all of the gene names which can be copy and pasted into the GeneMania cytoscape app, data frame of the names of the genes
#' @export
#'
#' @examples
#' cccn.cfn.tools:::ex.MakeDBInput(ex.nodenames)
MakeDBInput <- function(cccn_matrix, file.path.name = "db_nodes.txt") {
  utils::write.table(rownames(cccn_matrix), file = file.path.name, row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# Pulls nodenames from the cccn_matrix
#
# This helper function pulls the gene names from the cccn_matrix into a list 'nodenames'
#
# @param cccn_matrix dataframe of dataframes that represent the common clusters from the three distance calculations' clusters
# @return data frame of the names of the genes
cccn_to_nodenames <- function(cccn_matrix, nodenames.name = 'nodenames'){

  gene_names <- unique(rownames(cccn_matrix))

  nodenames <- data.frame(Gene.Names = gene_names, stringsAsFactors = FALSE)

  #return :)
  assign(nodenames.name, nodenames, envir = .GlobalEnv)
}


#' Get STRINGdb
#'
#' This function finds protein-protein interaction edges by consulting utilizing the STRINGdb database.
#'
#' @param cccn_matrix dataframe of dataframes that represent the common clusters from the three distance calculations' clusters
#' @param STRINGdb.name desired name for the output STRINGdb data frame; defaults to "string.edges"
#' @param nodes.name desired name for the output containing the list of all of the gene names; defaults to "nodenames"
#'
#' @return A data frame of edges from STRINGdb showing interactions between genes
#' @export
#'
#' @examples
#' cccn.cfn.tools:::ex.FindPPIEdges(ex.cccn_matrix, "ex.string.edges", "ex.nodenames")
GetSTRINGdb <- function(cccn_matrix, STRINGdb.name = "string.edges", nodenames.name = "nodenames") {
  cccn_to_nodenames(cccn_matrix, nodenames.name)

  if (!exists("STRINGdb")){                          # check if stringdb is libraried

    if(system.file(package="STRINGdb") == ""){       # check if stringdb is installed at all

      BiocManager::install("STRINGdb")               # install!
    }
    library(STRINGdb)                                # library it
  }

  # Initialize the STRING database object
  string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=0, input_directory="")

  # Retrieve the proteins from the STRING database
  string_proteins <- string_db$get_proteins()

  if (!"Gene.Names" %in% colnames(nodenames)) {
    stop("Column 'Gene.Names' not found in nodenames.")
  }

  # Map the genes to STRING IDs
  # please note that nodenames replaces the previous "input_dataset"; nodenames appears to work well :)
  # Gene.Names also replaces experimental
  mapped_genes <- string_db$map(nodenames, "Gene.Names", removeUnmappedRows = TRUE)
  print(utils::head(mapped_genes))

  # Retrieve the interactions for the mapped genes
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)

  # Convert protein IDs to gene names
  interactions$Gene.1 <- sapply(interactions$from, function(x) string_proteins[match(x, string_proteins$protein_external_id), "preferred_name"])
  interactions$Gene.2 <- sapply(interactions$to, function(x) string_proteins[match(x, string_proteins$protein_external_id), "preferred_name"])

  # Create the final edges dataframe from STRINGdb
  combined_edges <- interactions[, c("Gene.1", "Gene.2", "combined_score")]
  colnames(combined_edges) <- c("Gene.1", "Gene.2", "STRINGdb.combined_score")

  # assign
  assign(STRINGdb.name, combined_edges, envir = .GlobalEnv)
}


#' Process GM Edge File
#'
#' This function processes the GM edgefile and translates it back into gene names using the nodetable.
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

  split_names <- strsplit(edgetable$name, "\\|")                                                       # split the names up bc it has three pieces of info there
  edgetable$Gene.1 <- sapply(split_names, function(x)x[1])                                             # take the first thing; first ID name
  edgetable$Gene.1 <- sapply(edgetable$Gene.1, function(x)nodetable$query.term[nodetable$name == x])   # turn the ID into the gene name!
  edgetable$Gene.2 <- sapply(split_names, function(x)x[2])                                             # take the first thing; first ID name
  edgetable$Gene.2 <- sapply(edgetable$Gene.2, function(x)nodetable$query.term[nodetable$name == x])   # turn the ID into the gene name!

  edges <- edgetable[, c("Gene.1", "Gene.2", "normalized.max.weight")]      # sort into a new table with only our information and in the order we want

  colnames(edges) <- c("Gene.1", "Gene.2", "GM.weights")                    # rename so when put with everything, clearer where came from

  keep <- edges$Gene.1 %in% nodenames & edges$Gene.2 %in% nodenames         # which rows are we keeping
  all.edges <- edges[keep,]                                                 # copy 'em over

  assign(gm.network.name, as.data.frame(all.edges), envir = .GlobalEnv)                 # assign :)
}
