#' Find PPI Edges
#'
#' This function finds protein-protein interaction edges by combining STRINGdb and GeneMANIA databases.
#'
#' @param cccn_matrix dataframe of dataframes that represent the common clusters from the three distance calculations' clusters
#' @param db_filepaths vector of filepaths to data from other databases; defaults to empty vector
#' @param gm.network network produced by GeneMANIA if that was used; defaults to NA
#' @param ppi.network.name desired name of the output ppi network; defaults to "ppi.network"
#'
#' @return A data frame of combined edges from STRINGdb and provided database entries
#' @export
#'
#' @examples
#' gmpath <- "genemania-interactions.txt"
#' pack <- "cccn.cfn.tools"
#' gmfile <- system.file("genemania", gmpath, package = pack, mustWork = TRUE)
#' cccn.cfn.tools:::ex.FindPPIEdges(ex.cccn_matrix, ppi.network.name = "ex.ppi_network")
FindPPIEdges <- function(string.edges = NA, db_filepaths = c(), gm.network = NA, ppi.network.name = "ppi.network") {

  if (is.na(string.edges) && is.na(gm.network) && (length(db_filepaths) == 0)){
    stop("No data input.")
  }

  ppi.network <- data.frame()

  if(!is.na(string.edges)){
    if(!is.data.frame(string.edges)){
      stop("string.edges is not a dataframe; invalid formatting")
    }
    ppi.network <- rbind(ppi.network, string.edges)
  }

  if (!is.na(gm.network)){
    if(!is.data.frame(string.edges)){
      stop("gm.network is not a dataframe; invalid formatting")
    ppi.network <- rbind(ppi.network, gm.network)
    }
  }

  # Combine STRINGdb and GeneMANIA edges if gm_edges exists
  if(length(db_filepaths) != 0){
    for(path in db_filepaths){
      db_edges <- utils::read.table(path)
      ppi.network <- rbind(ppi.network, db_edges)
    }
  }

  assign(ppi.network.name, ppi.network, envir = .GlobalEnv)

}
