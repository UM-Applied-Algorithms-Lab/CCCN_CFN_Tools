#' Creates sample data that mimics the format of global variables eu_ptms_list, sp_ptms_list, and sed_ptms_list
#' 
#' This function adds sample data: Fruits1, Fruits2, and Fruits3 to the global enviroment, edit as desired. 
#' @param save Will save the lists as files in the same directory as SampleDataCreator.r
#'
#' @examples
#' Fruitify()
#' Fruitify(TRUE)
Fruitify <- function(save = FALSE){
  
#Database
#Use this space to keep track of the fruits in the data
###
#' "Apple", "Avocado", "Banana", "Blackberry", "Blueberry", "Breadfruit", "Catmon",
#' "Cherry", "Citron", "Date", "Dragonfruit", "Fig", "Guava", "Kiwi", "Kumquat", "Lemon", "Lime",
#' "Lychee", "Mango", "Melon", "Orange", "Papaya", "Peach", "Pear", "Persimmon", "Pineapple",
#' "Plum", "Pomegranate", "Raspberry", "Starfruit", "Strawberry", "Tayberry"
###

  
  
#EDIT DATA HERE
#Fruits 1 data - Add, remove, change as needed

Fruits1.Data <- list(
  c("Apple", "Avocado", "Banana", "Blackberry", "Blueberry"),
  c("Breadfruit", "Catmon", "Cherry", "Citron", "Date", "Dragonfruit", "Fig"),
  c("Guava", "Kiwi", "Kumquat", "Lemon"),
  c("Lychee", "Mango", "Melon", "Orange", "Papaya", "Peach", "Pear", "Persimmon", "Pineapple", "Plum", "Pomegranate", "Raspberry", "Starfruit", "Strawberry", "Tayberry")
)

#Fruits 2 data - Add, remove, change as needed 
Fruits2.Data <- list(
  c("Blueberry", "Breadfruit", "Catmon"),
  c("Cherry", "Citron", "Date"),
  c("Apple", "Avocado", "Banana", "Blackberry"),
  c("Dragonfruit", "Fig", "Guava"),
  c("Kiwi", "Kumquat", "Lemon"),
  c("Lime", "Lychee", "Mango"),
  c("Orange", "Papaya", "Peach"),
  c("Persimmon", "Pineapple", "Plum"),
  c("Pomegranate", "Raspberry"),
  c("Starfruit", "Strawberry", "Tayberry")
)

#Fruits 3 data - Add, remove, change as needed
Fruits3.Data <- list(
  c("Cherry", "Citron", "Date", "Dragonfruit", "Fig"),
  c("Blueberry", "Breadfruit", "Catmon"),
  c("Guava", "Kiwi", "Kumquat", "Persimmon", "Pineapple", "Plum", "Pomegranate", "Raspberry", "Starfruit", "Strawberry", "Tayberry"),
  c("Apple", "Avocado", "Banana", "Blackberry", "Lemon", "Lime"),
  c("Lychee", "Mango", "Melon", "Orange", "Papaya", "Peach", "Pear")
)

#If you want to add more data, add it here
#Fruits4.Data <- list(vector, vector, vector)



#Anything you added must be put into this list (like Fruits4.Data)! Run Fruitify() after sourcing this file to get the lists added to namespace
Total <- list(Fruits1.Data, Fruits2.Data, Fruits3.Data)



#DO NOT EDIT - Processing code
for(a in 1:length(Total)){
  count <- 0L
  for(b in 1:length(Total[[a]])){
    Total[[a]][[b]] <- as.data.frame(Total[[a]][[b]])
    colnames(Total[[a]][[b]])<- c("Gene.Name")
    count <- count + 1L
    Total[[a]][[b]]$groups <- rep(count, length(Total[[a]][[b]]))
  }
  string <- paste("Fruits", as.character(a), "_ptm_list", sep="") #Change these names to be more descriptive
  temp <- Total[[a]]
  assign(string, temp, envir = .GlobalEnv)
  if(save) save(temp, file = paste(string, ".rda", sep="")) #Save the data if true - BROKEN
}

}