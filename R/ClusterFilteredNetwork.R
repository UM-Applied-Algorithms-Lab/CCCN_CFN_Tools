
# Function to extract gene names from peptide names
pepgene <- function(peps) {
  unique(sapply(peps, function(x) unlist(strsplit(x, " ", fixed=TRUE))[1])) #Definitely a cry for help. 100%. You can't tell me it isn't. 6 NESTED FUNCTIONS.  
}