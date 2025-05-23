#Setup
if(!exists("eu_ptms_list")){      #Check if global variables are already made as to not run MakeClusterList multiple times
  set.seed(1)                     #Set the seed (very important)
  
  #Load Sample data
  suppressWarnings(
    if(!exists("ptmtable")){
      try(load("CCCN_CFN_Tools/data/ptmtable.rda"), silent=TRUE)
      try(load("data/ptmtable.rda"), silent=TRUE)
      try(load("../data/ptmtable.rda"), silent=TRUE)
      try(load("../../data/ptmtable.rda"), silent=TRUE)
      if(!exists("ptmtable")){stop("Cannot find ptmtable in CCCN_CFN_Tools, please be in the CCCN_CFN_Tools directory and make sure ptmtable exists!")}
    })
  

  #Make Global Variables
  sink("noprint")                 #Suppress print statements from function
  MakeClusterList(ptmtable)       #Create sample data - #BUG - writes 'species scores not available' (dont worry about this for now)
  sink()
}

#Unit Testing - Test cluster sizes MAKE THIS A BETTER TEST (if any, this doesn't seem to do anything unique, I'd like to test to see if it writing plots or not)
test_that("Testing eu_ptms_list", {expect_equal(length(sapply(eu_ptms_list, function(x) dim(x)[1])), 6)})
test_that("Testing sed_ptms_list", {expect_equal(length(sapply(sed_ptms_list, function(x) dim(x)[1])), 6)})
test_that("Testing sp_ptms_list", {expect_equal(length(sapply(sp_ptms_list, function(x) dim(x)[1])), 88)})


#Clean Up
if(file.exists("noprint")) file.remove("noprint") #Clean up file created by sink()
