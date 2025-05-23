#' Creates sample data that mimics the format of global variables eu_ptms_list, sp_ptms_list, and sed_ptms_list
#' 
#' This function adds sample data: Fruits1, Fruits2, and Fruits3 to the global enviroment, edit as desired. 
#' @param save Will save the lists as files in the same directory as SampleDataCreator.r
#'
#' @examples
#' Fruitify()
#' Fruitify(TRUE)
Fruitify <- function(save = FALSE){
  
  #DO NOT EDIT 
  counter <- 0 
  makegroup <- function(Gene.Name){counter <- counter+1; return(rep(counter, length(Gene.Name)))} 
  #DO NOT EDIT

#Create some sample data
Fruits1 <- list(
    data.frame(
      Gene.Name = c("Apple", "Orange", "Banana", "Tomato", "Blueberry"), #Edit me
      group = rep(1L, 5)
    ),
    
    data.frame(
      Gene.Name = c("Apricot", "Avocado", "Blackberry", "Mango", "Eggplant"), #Edit me
      group = rep(2L, 5)
    ),
    
    data.frame(
      Gene.Name = c("Cherry", "Cucumber", "Strawberry", "Blueberry", "Green Apple"), #Edit me
      group = rep(3L, 5)
    )
  
)

#Create more sample data
Fruits2 <- list(
  data.frame( 
      Gene.Name = c("Cherry", "Strawberry"),
      group = rep(1L, 2)
    ),
  
  data.frame(
    Gene.Name = c("Blueberry", "Green Apple"),
    group = rep(2L, 2)
  ),
  
  data.frame(
    Gene.Name = c("Apple", "Orange", "Banana", "Tomato", "Apricot", "Avocado", "Blackberry"),
    group = rep(3L, 7)
  ),
  
  data.frame(
    Gene.Name = c("Mango", "Eggplant"),
    group = rep(4L, 2)
  )
)

#Create more sample data
Fruits3 <- list(
  data.frame(
    Gene.Name = c("Apple", "Tomato", "Orange", "Avocado", "Blueberry", "Banana"),
    group = rep(1L, 6)
  ),
  
  data.frame(
    Gene.Name = c("Blackberry", "Cherry", "Cucumber"),
    group = rep(2L, 3)
  ),
  
  data.frame(
    Gene.Name = c("Apricot", "Mango", "Eggplant"),
    group = rep(3L, 3)
  ),
  
  data.frame(
    Gene.Name = c("Strawberry", "Blueberry", "Green Apple"),
    group = rep(4L, 3)
  )
)

assign("fruits", Fruits1, envir = .GlobalEnv)   #Assigns sample data to global variable 
assign("fruits2", Fruits2, envir = .GlobalEnv)
assign("fruits3", Fruits3, envir = .GlobalEnv)  

if(save) save(Fruits1, file = "Fruits1.rda")   #Save the data to file if TRUE was passed in
if(save) save(Fruits2, file = "Fruits2.rda")
if(save) save(Fruits2, file = "Fruits3.rda")
}
