
Fruitify <- function(save = FALSE){

#Create some sample data
Fruits <- list(
    data.frame(
      Fruits.Name = c("Apple", "Orange", "Banana", "Tomato", "Blueberry"),
      group = rep(1L, 5)
    ),
    
    data.frame(
      Fruits.Name = c("Apricot", "Avocado", "Blackberry", "Mango", "Eggplant"),
      group = rep(2L, 5)
    )
  
)

#Create more sample data
Fruits2 <- list(
  data.frame( 
      Fruits.Name = c("Apple", "Orange"),
      group = rep(1L, 2)
    ),
  
  data.frame(
    Fruits.Name = c("Banana", "Tomato"),
    group = rep(2L, 2)
  ),
  
  data.frame(
    Fruits.Name = c("Apricot", "Avocado", "Blackberry"),
    group = rep(3L, 3)
  ),
  
  data.frame(
    Fruits.Name = c("Mango", "Eggplant"),
    group = rep(4L, 2)
  )
)

#Create more sample data
Fruits3 <- list(
  data.frame(
    Fruits.Name = c("Apple", "Avocado", "Blueberry", "Banana"),
    group = rep(1L, 4)
  ),
  
  data.frame(
    Fruits.Name = c("Blackberry", "Tomato", "Orange"),
    group = rep(2L, 3)
  ),
  
  data.frame(
    Fruits.Name = c("Apricot", "Mango", "Eggplant"),
    group = rep(3L, 3)
  )
)

assign("fruits", Fruits, envir = .GlobalEnv)   #Assigns sample data to global variable 
assign("fruits2", Fruits2, envir = .GlobalEnv)
assign("fruits3", Fruits3, envir = .GlobalEnv)  

if(save) save(Fruits, file = "listtest.rda")   #Save the data to file if TRUE was passed in
if(save) save(Fruits2, file = "listtest2.rda")
if(save) save(Fruits2, file = "listtest3.rda")
}
