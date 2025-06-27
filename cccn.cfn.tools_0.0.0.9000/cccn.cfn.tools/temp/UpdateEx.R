#' Updates (most) data in the data directory
#' 
#' @param data Use if you don't want to use the standard example ptmtable
#' @examples
#' UpdateEx()
UpdateEx <- function(table = ex.ptmtable){
  if(!"verify12345.txt" %in% list.files("temp")) stop("Verification file not found. Please check your directory.")
  setwd("data") #Temporarily save to data
  set.seed(1)
  
  #~MakeClusterList example data update~
  MakeClusterList(table) #Run
  
  #Name data
  ex.ptm.correlation.matrix <- ptm.correlation.matrix
  ex.eu_ptms_list <- eu_ptms_list
  ex.sp_ptms_list <- sp_ptms_list
  ex.sed_ptms_list <- sed_ptms_list
  
  #Save data to files
  save(ex.eu_ptms_list, file="eu_ptms_list.rda")
  save(ex.sp_ptms_list, file="sp_ptms_list.rda")
  save(ex.sed_ptms_list, file="sed_ptms_list.rda")
  save(ex.ptm.correlation.matrix, file="ptm_correlation_matrix.rda")
  
  
  #~MakeCorrelationNetwork example data update~
  MakeCorrelationNetwork(1) #Run
  
  #Name data
  ex.cccn_matrix <- cccn_matrix
  
  #Save data to files
  save(ex.cccn_matrix, file="cccn_matrix.rda")
  
  
  #PPI edges
  #to do
  
  
  #ClusterFilteredNetwork
  #to do
  
  
  
  setwd("..") #Go back to CCCN_Tools
  
}