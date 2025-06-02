# Note for Maddie bc Ik you're gonna forget
# STRINGdb works, cccn_to_nodenames works (make sure cccn is a dataframe I think; igraph object seems silly to work with)
# get.GM.edges is the next thing to be fixed; look at making GeneMania work on the command line
# after that, this file should be done


#' Loads and filters the GeneMania file given a vector of gene names.Add commentMore actions
#'
#' This helper function loads the GeneMania file
#' and filters out required nodes.
#'
#' @param nodenames A vector containing the names of the relevant genes.
#' @param gmfilename The path to the GeneMania file initialized to "genemania-interactions.txt".
#' @return A data frame with the relevant GeneMania data.
#' @export
#'
#' @examples
#' get.GM.edgefile(nodenames, gmfilename)
get.GM.edgefile <- function(nodenames, gmfilename = "genemania-interactions.txt"){

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
  #Removes the column "Network" that just tells what paper this was published in and assigns it to the namespace
  assign('gm_edges', gm_edges[ , -5])
}

#' Pulls nodenames from the cccn_matrix
#'
#' This helper function pulls the gene names from the cccn_matrix into a list 'nodenames'
#'
#' @param cccn_matrix dataframe of dataframes that represent the common clusters from the three distance calculations' clusters
#' @return data frame of the names of the genes
#' @export
#'
#' @examples
#' cccn_to_nodenames(cccn_matrix)
cccn_to_nodenames <- function(cccn_matrix){
  # initialize as an empty list
  nodenames <- data.frame(matrix(ncol = 1, nrow = 0))
  colnames(nodenames) <- "Gene.Names"
  len <- nrow(nodenames)

  # steal the row names
  cccn_rows <- rownames(cccn_matrix)

  #loop through row names
  for (i in 1:length(cccn_rows)){

    #assign the var genename to the name of the gene
    genename <- cccn_rows[i]

    #check if the genename is already in the list
    if (!(genename %in% nodenames)){

      #go to next entry (which doesn't exist yet)
      len <- len +1

      #create the entry
      nodenames[len,] <- genename

    }
  }

  #return :)
  assign("nodenames", nodenames)
}



#eu.sp.sed.ptms.data into something cccn_matrix - like:
#  for(i in 1:nrow(nodenames)){
#    name <- strsplit(nodenames[i, 1], ' ')[[1]][1]
#    nodenames[i,1] <- name
#  }



#' Make file for GeneMania input
#'
#' This function outputs a file the user can take to create the GeneMania edgefile.
#'
#' @param cccn_matrix dataframe of dataframes that represent the common clusters from the three distance calculations' clusters
#'
#' @return A file with all of the gene names which can be copy and pasted into the GeneMania website, data frame of the names of the genes
#' @export
#'
#' @examples
#' find_ppi_edges(cccn_matrix)
make_gm_input <- function(cccn_matrix) {
  cccn_to_nodenames(cccn_matrix)
  write.table(nodenames$Gene.Names, file = "gm_nodes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#' Find PPI Edges
#'
#' This function finds protein-protein interaction edges by combining STRINGdb and GeneMANIA databases.
#'
#' @param cccn_matrix dataframe of dataframes that represent the common clusters from the three distance calculations' clusters
#' @param gmfilename The path to the GeneMania file initialized to "genemania-interactions.txt".
#'
#' @return A data frame of combined edges from STRINGdb and GeneMANIA.
#' @export
#'
#' @examples
#' find_ppi_edges(cccn_matrix, gmfilename = 'genemania/genemania-interactions.txt')
find_ppi_edges <- function(cccn_matrix, gmfilename = "genemania-interactions.txt") {

  #test if nodenames already exists and therefore if they ran optional step 3
  if(!exists("nodenames")){
    #find the nodenames
    cccn_to_nodenames(cccn_matrix)} else{ #if nodenames does exist then they should have ran step 3 and have gm data
    # Get GeneMANIA edges (Reminder these ^ V are both assigned to the global)
    get.GM.edgefile(nodenames, gmfilename)
  }

  # Initialize the STRING database object
  string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=0, link_data="detailed", input_directory="")

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

  # Filter interactions based on evidence types
  #These were picked based on my limited knowledge of ppis
  #may easily be changed
  str.e <- interactions[interactions$experimental > 0, ]
  str.n <- interactions[interactions$neighborhood > 0, ]
  str.d <- interactions[interactions$database > 0, ]
  str.co <- interactions[interactions$cooccurence > 0, ]
  str.ce <- interactions[interactions$coexpression > 0, ]

  # Combine filtered interactions
  combined_interactions <- unique(rbind(str.e, str.n, str.d, str.co, str.ce))

  # Assign edge types
  combined_interactions$edgeType <- "STRINGdb"
  combined_interactions[combined_interactions$database > 0, "edgeType"] <- "database"
  combined_interactions[combined_interactions$neighborhood > 0, "edgeType"] <- "neighborhood"
  combined_interactions[combined_interactions$experimenal > 0, "edgeType"] <- "experimental"
  combined_interactions[combined_interactions$cooccurence > 0, "edgeType"] <- "cooccurence"
  combined_interactions[combined_interactions$coexpression > 0, "edgeType"] <- "coexpression"

  # Calculate weights
  combined_interactions$Weight <- rowSums(combined_interactions[, c("experimental", "database", "neighborhood", "cooccurence", "coexpression")])
  combined_interactions$Weight <- combined_interactions$Weight / 1000

  # Create the final edges dataframe from STRINGdb
  combined_edges <- combined_interactions[, c("Gene.1", "Gene.2", "Weight", "edgeType")]

  # Combine STRINGdb and GeneMANIA edges if gm_edges exists
  if(exists("gm_edges")){
    combined_ppi_network <- rbind(combined_edges, gm_edges)
    assign("ppi_network", combined_ppi_network)} else{ #if gm_edges does not exist then do not combine and only use those from STRINGdb
    assign("ppi_network", combined_edges)
  }
}
