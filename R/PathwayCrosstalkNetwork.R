#'Helper function to take information from a bioplanet file and format it in the desired way. Assumes that the function was downloaded from excel. 
#'
#'@param filename The name of the bioplanet pathway .csv file    
#'@return A formatted version of the bioplanet csv file
#'@export
#'
#'@examples
#'print("TO DO")
bioplanet.innit <- function(fn){
  print("Unused; may delete later")
}

#'
#'
#'@param filename The name of the bioplanet pathway .csv file    
#'@return No clue (The result is a matrix with pathway names in columns and individual clusters as rows.)
#'@export
#'
#'@examples
#'print("TO DO")
PathwayCrosstalkNetwork <- function(filename){
  if(!file.exists(filename)) stop(paste(filename, "not found. Plese check your working directory."))
  
  #Read file in, converts to dataframe like with rows like: PATHWAY_ID | PATHWAY_NAME | GENE_ID | GENE_SYMBOL
  bioplanet <- read.csv(filename, stringsAsFactors = F) 
  
  ###
  #Add excel fix here? 
  ###
  
  #Return
  return(bioplanet) 
  
}