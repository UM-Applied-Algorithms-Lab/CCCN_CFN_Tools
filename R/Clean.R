fix.excel <- function(vector) {  
  
  #hashtable: input what needs to be fixed, output the fixed version
  correct.to <- c("CDK1", "SEPT1", "SEPT2", "SEPT3", "SEPT4", "SEPT5", "SEPT7", "SEPT8", "SEPT9", "SEPT10", "SEPT11", "SEPT15", "SEPT6", "POU2F1", "POU2F2", "POU5F1", "POU5F1", "POU3F1", "POU3F2", "POU2F3", "MARCH1", "MARCH2", "MARCH3", "MARCH4", "MARCH5", "MARCH6", "MARCH7", "MARCH8", "MARCH9", "MARCH10", "MARCH11", "SMAP", "EFCAB13", "BRK1", "NYAP1", "LAMTOR1", 'AP1AR', "DEC1", "SEPT14")
  names(correct.to) <- c("CDC2", "1-Sep", "2-Sep", "3-Sep", "4-Sep", "5-Sep", "7-Sep", "8-Sep", "9-Sep", "10-Sep", "11-Sep", "15-Sep", "6-Sep", "1-Oct", "2-Oct", "3-Oct", "4-Oct", "6-Oct", "7-Oct", "11-Oct", "1-Mar", "2-Mar", "3-Mar", "4-Mar", "5-Mar", "6-Mar", "7-Mar", "8-Mar", "9-Mar", "10-Mar", "11-Mar", "C11orf58", 'C17orf57', 'C3orf10',  'C7orf51', "C11orf59", "C4orf16", "1-Dec", "14-Sep")
  
  indeces <- vector %in% fixme #indecies that need to be fixed
  temp.vec <- vector[indeces] #create a temp vec
  
  for(i in 1:length(temp.vec)) temp.vec[[i]] <- correct.to[[temp.vec[[i]]]] #Correct it
  
  vector[indeces] <- temp.vec #Replace everything with temp vec
  
  return(vector)
}