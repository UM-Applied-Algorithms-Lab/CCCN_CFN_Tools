#' Converts Bioplanet pathways from (<https://tripod.nih.gov/bioplanet/>)  into a list of pathways whose elements are the genes in each pathway. Edge weights are either the PTM Cluster Weight or according to the Jaccard Similarity.
#' 
#' @param filename The name of the bioplanet pathway .csv file    
#' @return No clue (The result is a matrix with pathway names in columns and individual clusters as rows.)
#' @export
#'
#' @examples
#' print("TO DO")
PathwayCrosstalkNetwork <- function(filename, PCNname = "pcn_matrix"){
  #Error catch if file/directory is wrong
  if(!file.exists(filename)) stop(paste(filename, "not found. Plese check your working directory."))
  
  #Read file in, converts to dataframe like with rows like: PATHWAY_ID | PATHWAY_NAME | GENE_ID | GENE_SYMBOL
  bioplanet <- read.csv(filename, stringsAsFactors = F) 
  
  
}