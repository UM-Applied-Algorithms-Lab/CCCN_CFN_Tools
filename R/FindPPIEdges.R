#' Find PPI Edges
#'
#' This function finds protein-protein interaction edges by combining STRINGdb and GeneMANIA databases.
#'
#' @param cccn.matrix dataframe of dataframes that represent the common clusters from the three distance calculations' clusters
#' @param db.filepaths vector of filepaths to data from other databases; defaults to empty vector
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
#' cccn.cfn.tools:::ex.FindPPIEdges(ex.cccn.matrix, ppi.network.name = "ex.ppi.network")
FindPPIEdges <- function(string.edges = NA, gm.network = NA, db.filepaths = c(), ppi.network.name = "ppi.network") {

  if (!is.data.frame(string.edges) && !is.data.frame(gm.network) && (length(db.filepaths) == 0)){
    stop("No data input.")
  }

  ppi.network <- data.frame()

  bind_ppis <- function(ppi.network, edgefile){

    columname <- colnames(edgefile)[3]            # name the column
    ppi.network[[columname]] <- NA                # make new column
    len <- length(colnames(ppi.network))          # record num cols

    for(i in 1:length(rownames(edgefile))){       # iterate through rows

      gene1 <- edgefile[i,1]
      gene2 <- edgefile[i,2]
      weight <- edgefile[i,3]

      row <- which(ppi.network[[1]] == gene1 & ppi.network[[2]] == gene2)

      if(length(row) == 0){
        row <- which(ppi.network[[1]] == gene2 & ppi.network[[2]] == gene1)

        if(length(rownames(row)) == 0){                             # if this combo DNE in ppi

          newrow <- as.data.frame(matrix(NA, nrow = 1, ncol = len)) # make a new dataframe with our info
          colnames(newrow) <- colnames(ppi.network)                 # make the column names the same
          newrow[[1]] <- gene1                                      # add the data
          newrow[[2]] <- gene2
          newrow[[len]] <- weight

          ppi.network <- rbind(ppi.network, newrow)                 # combine

          next # DON'T DOUBLE ASSIGN, DUMMY

                  }
      }

      for(index in row){
      ppi.network[index, len] <- weight
      }

    }

    return(ppi.network)

    }

  # ppi is essentially initialized to the string.edges
  # TODO: FIX FOR CASES WHERE USER DOES NOT USE STRINGDB!!!!
  if(is.data.frame(string.edges)){
    ppi.network <- rbind(ppi.network, string.edges)
  }

  # Combine STRINGdb and GeneMANIA edges if gm.network exists
  if (is.data.frame(gm.network)){
    ppi.network <- bind_ppis(ppi.network, gm.network)
  }

  # bind all of the db edges to the ppi.network
  if(length(db.filepaths) != 0){
    for(path in db.filepaths){
      db.edges <- utils::read.table(path)
      ppi.network <- bind_ppis(ppi.network, db.edges)
    }
  }

  assign(ppi.network.name, ppi.network, envir = .GlobalEnv)

}
