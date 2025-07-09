#' Build PPI Network
#'
#' Uses all inputted database information on protein-protein interactions in order to create a combined PPI-network.
#'
#' @param stringdb.edges Data frame of consisting of the network of interactions from the genes of study pulled from the STRINGdb database
#' @param gm.network GeneMANIA network of protein-protein interactions from the genes of study; defaults to NA
#' @param db.filepaths A vector of paths to the additional ppi network files; defaults to an empty vector
#' @param ppi.network.name Desired name for the output protein-protein interaction network using all entered database input; defaults to ppi.network
#'
#' @return A dataframe representing how strongly proteins are known to interact based on previous research contained within provided database networks
#' @export
#'
#' @examples
#' ex.string <- ex.stringdb.edges
#' ex.gm <- ex.gm.network
#' BuildPPINetwork(stringdb.edges = ex.string, gm.network = ex.gm, ppi.network.name = "ex.ppi.network")
#' utils::head(ex.ppi.network)
BuildPPINetwork <- function(stringdb.edges = NA, gm.network = NA, db.filepaths = c(), ppi.network.name = "ppi.network") {

  if (!is.data.frame(stringdb.edges) && !is.data.frame(gm.network) && (length(db.filepaths) == 0)){
    stop("No data input.")
  }

  ppi.network <- data.frame()

  bind_ppis <- function(ppi.network, edgefile){

    columname <- colnames(edgefile)[4]            # name the column
    ppi.network[[columname]] <- NA                # make new column
    len <- length(colnames(ppi.network))          # record num cols

    for(i in 1:length(rownames(edgefile))){       # iterate through rows

      gene1 <- edgefile$Gene.1[i]
      gene2 <- edgefile$Gene.2[i]
      int.type <- edgefile[i,3]
      weight <- edgefile[i,4]

      row <- which(ppi.network[[1]] == gene1 & ppi.network[[2]] == gene2)

      if(length(row) == 0){
        row <- which(ppi.network[[1]] == gene2 & ppi.network[[2]] == gene1)
      }

      if(length(rownames(row)) == 0){                               # if this combo DNE in ppi

          newrow <- as.data.frame(matrix(NA, nrow = 1, ncol = len)) # make a new dataframe with our info
          colnames(newrow) <- colnames(ppi.network)                 # make the column names the same
          newrow[[1]] <- gene1                                      # add the data
          newrow[[2]] <- gene2
          newrow[[3]] <- int.type
          newrow[[len]] <- weight

          ppi.network <- rbind(ppi.network, newrow)                 # combine

        } else {                                  # if combo DOES exist

        for(index in row){                      # iterate through the combinations
        ppi.network[index, len] <- weight       # assign the weight

        ppi.network$Interaction.Type[index] <- paste(ppi.network$Interaction.Type[index], int.type, sep = ", ")

        }
      }
    }

    maxi <- max(ppi.network[[columname]], na.rm = TRUE)                                   # get max of column
    ppi.network[[columname]] <- sapply(ppi.network[[columname]], function(x) (x / maxi))  # divide by max to get scale of 0-1

    return(ppi.network)

    }

  # ppi is essentially initialized to the stringdb.edges if stringdb.edges isn't NA
  if(is.data.frame(stringdb.edges)){
    if("Gene.1" %in% colnames(stringdb.edges)[-3] & "Gene.2" %in% colnames(stringdb.edges)[-3] & length(colnames(stringdb.edges)) == 4){  # check formatting

      columname <- colnames(stringdb.edges)[3]                                              # get column name

      ppi.network <- rbind(ppi.network, stringdb.edges)                                     # initialize ppi.network w rbind

      maxi <- max(ppi.network[[columname]], na.rm = TRUE)                                   # get the max
      ppi.network[[columname]] <- sapply(ppi.network[[columname]], function(x) (x / maxi))  # divide each by max for scale 0-1
    }
    } else {
    print("Improper formatting of stringdb.edges. Ensure creation by GetSTRINGdb.")         # warning message else
  }

  # Combine STRINGdb and GeneMANIA edges if gm.network exists
  if (is.data.frame(gm.network)){
    if("Gene.1" %in% colnames(gm.network)[-3] & "Gene.2" %in% colnames(gm.network)[-3] & length(colnames(gm.network)) == 4){  # check formatting
      if (length(rownames(ppi.network)) == 0){                                                                                # initialize if good
        ppi.network <- rbind(ppi.network, gm.network)                                                                         # and no STRING
      } else{                                                                                                                 # bind w my func
      ppi.network <- bind_ppis(ppi.network, gm.network)                                                                       # if STRING yes
      }
    } else {                                                                                                                  # but if format BAD
      print("Improper formatting of gm.network. Ensure processing by ProcessGMEdgefile.")                                     # warning message
    }
  }

  # bind all of the db edges to the ppi.network
  if(length(db.filepaths) != 0){                         # if they entered database files of their own
    for(path in db.filepaths){                           # iterate through the file paths

      db.edges <- utils::read.table(path)                                                                          # get the info from the file
      if(!("Gene.1" %in% colnames(db.edges)[-3]) | !("Gene.2" %in% colnames(db.edges)[c(-3, -4)])){                       # check if colnames correct
        print("Improper naming of column names in the following file:")                                            # warning message if not
        cat(path)
        print("First two columns should be labeled 'Gene.1' and 'Gene.2'. Skipping this database")
        next                                                                                                       # skip if improper colnames
      }
      if(!(length(colnames(db.edges)) == 4)){                                                                      # check if wrong num cols
        print("Improper number of columns in the following file:")                                                 # warning message if not
        cat(path)
        print("There should be four columns with the third being the interaction tybe and the fourth being the weights of the edges")
        next                                                                                                       # skip if num cols
      }
      if (length(rownames(ppi.network)) == 0){           # if didn't use STRINGdb or GeneMANIA it's a simple Rbind
        ppi.network <- rbind(ppi.network, db.edges)
      } else{                                            # did use STRINGdb and/or GeneMANIA it's my binding function
        ppi.network <- bind_ppis(ppi.network, db.edges)
      }
    }
  }

  ppi.network <- unique(ppi.network)
  assign(ppi.network.name, ppi.network, envir = .GlobalEnv)
}
