#draft of list.common

#TODO: Find some better way of appending, improve time complexity?, Fruits to Gene issue - data must have a column named Gene.Names, workaround?

# Helper function to find intersections of clusters
#'
#' Finds common elements between clusters in three lists.
#'
#' @param list1 A list of ptms.
#' @param list2 A list of ptms.
#' @param list2 A list of ptms.
#' @param keeplength Minimum size of intersections to keep.
#' @return A list of common clusters.
#' @examples
#' list.common(list1, list2, list3 keeplength = 3)
draft <- function(list1, list2, list3, keeplength = 3, colname = "Gene.Name"){
  
  #Make sure that the desired column exists in sublists of lists
  if(!TRUE %in% sapply(list1, function(x) "Gene.Name" %in% names(x))) stop("List1 does not have a Gene.Name column")
  if(!TRUE %in% sapply(list2, function(x) "Gene.Name" %in% names(x))) stop("List2 does not have a Gene.Name column")
  if(!TRUE %in% sapply(list3, function(x) "Gene.Name" %in% names(x))) stop("List3 does not have a Gene.Name column")
  
  #Convert lists into groups of ptms
  list1.ptms <- lapply(list1, function(x){x$"Gene.Name"}) #These are lists of character vectors
  list2.ptms <- lapply(list2, function(y){y$"Gene.Name"}) #CHANGE "Fruits" to "Gene" or figure out some workaround
  list3.ptms <- lapply(list3, function(z){z$"Gene.Name"})
  
  #Find all the matching intersections of list1 and list2
  returnme <- list()  #Create an empty list to hold those intersections
  
  for(a in 1:length(list1.ptms)){
    for(b in 1:length(list2.ptms)){
      for(c in 1:length(list3.ptms)){
        temp <- Reduce(intersect, list(list1.ptms[[a]], list2.ptms[[b]], list3.ptms[[c]])) #Take the intersection of 3 character vectors
        if(length(temp) > keeplength) returnme[[length(returnme)+1]] <- temp               #And only add it if it has enough values
      }
    }
  }
  return(returnme) 
}


#Code Graveyard - REMOVE ME 
#finalsect <- list()
#Find all the matching intersections between tempsect and list3.ptms
#for(j in 1:length(tempsect)){
#  for(k in 1:length(list3.ptms)){
#    temp <- intersect(tempsect[[j]], list3.ptms[[k]])
#    if(length(temp) > keeplength) finalsect[[length(finalsect)+1]] <- temps
#  }
#}

#temp <- intersect(list1.ptms[[a]], list2.ptms[[b]])  #Find the intersections of character vector A and all character vectors in list2.ptms
#if(length(temp) > keeplength) tempsect[[length(tempsect)+1]] <- temp #CHANGE THIS HORRIBLE MONSTROSITY ASAP