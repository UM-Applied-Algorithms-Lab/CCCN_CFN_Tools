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

#' Find PPI Edges
#'
#' This function finds protein-protein interaction edges by combining STRINGdb and GeneMANIA databases.
#'
#' @param cccn_matrix dataframe of dataframes that represent the common clusters from the three distance calculations' clusters
#' @param db_filepaths vector of filepaths to data from other databases; defaults to empty vector
#'
#' @return A data frame of combined edges from STRINGdb and provided database entries
#' @export
#'
#' @examples
#' gmfile <- system.file("genemania", "genemania-interactions.txt", package = "cccn.cfn.tools", mustWork = TRUE)
#' cccn.cfn.tools:::ex.find_ppi_edges(ex.cccn_matrix)
FindPPIEdges <- function(cccn_matrix, db_filepaths = c(), gm.network = NA, ppi.network.name = "ppi.network") {

  cccn_to_nodenames(cccn_matrix)

  # Initialize the STRING database object
  string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=0, input_directory="")

  # Retrieve the proteins from the STRING database
  string_proteins <- string_db$get_proteins()
  print(dim(string_proteins))

  if (!"Gene.Names" %in% colnames(nodenames)) {
    stop("Column 'Gene.Names' not found in nodenames.")
  }

  # Map the genes to STRING IDs
  #please note that nodenames replaces the previous "input_dataset"; nodenames appears to work well :)
  #Gene.Names also replaces experimental
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

  if (!is.na(gm.network)){
    combined_edges <- rbind(combined_edges, gm.network)
  }

  # Combine STRINGdb and GeneMANIA edges if gm_edges exists
  if(length(db_filepaths) != 0){
    combined_ppi_network <- combined_edges
    for(path in db_filepaths){
      db_edges <- utils::read.table(path)
      combined_ppi_network <- rbind(combined_ppi_network, db_edges)
    }
    assign(ppi.network.name, combined_ppi_network, envir = .GlobalEnv)} else{ #if db_edges does not exist then do not combine and only use those from STRINGdb
    assign(ppi.network.name, combined_edges, envir = .GlobalEnv)
  }
}
