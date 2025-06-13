#' Converts Bioplanet pathways from (<https://tripod.nih.gov/bioplanet/>)  into a list of pathways whose elements are the genes in each pathway. Edge weights are either the PTM Cluster Weight or according to the Jaccard Similarity.
#' 
#' @param file The name of the bioplanet pathway .csv file    
#' @return No clue (The result is a matrix with pathway names in columns and individual clusters as rows.)
#' @export
#'
#' @examples
#' print("TO DO")
PathwayCrosstalkNetwork <- function(file, PCNname = "pcn_matrix"){
  #Error catch if file/directory is wrong
  if(!file.exists(file)) stop(paste(file, "not found. Plese check your working directory."))
  
  #Read file in, converts to dataframe like with rows like: PATHWAY_ID | PATHWAY_NAME | GENE_ID | GENE_SYMBOL
  bioplanet <- read.csv(file, stringsAsFactors = F) 
  
  #Turn bioplanet into a list of pathways. Pathways are character vectors comprised of gene names 
  pathways.list <- plyr::dlply(bioplanet, plyr::.(PATHWAY_NAME)) #Turn the bioplanet .csv into a list of data frames. Each data frame stores genes with the same PATHWAY_ID 
  pathways.list <- lapply(bioplanet, `[`, -c(1:3)) #Removes PATHWAY_ID | PATHWAY_NAME | GENE_ID from every single data frame in the list
  pathways.list <- lapply(pathways.list, unlist, use.names = FALSE) #Since data frames are 1 row, turn data frames into character vectors
  
  
}