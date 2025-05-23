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

  sink("noprint")                 #Suppress print statements from function
  MakeClusterList(ptmtable)       #Create sample data - #BUG - writes 'species scores not available' (dont worry about this for now)
  sink()
}

#Unit Tests for the global variables (lists of dataframes) created by MakeClusterList. Check "Gene" Names.
#eu_ptms_list
test_that("Testing eu_ptms_list Cluster 1, Gene 2", {expect_equal(eu_ptms_list[[1]]$Gene.Name[2], "ABCA1 ubi K2023")})
test_that("Testing eu_ptms_list Cluster 3, Gene 4", {expect_equal(eu_ptms_list[[3]]$Gene.Name[4], "ABCB6 ubi K482")})
test_that("Testing eu_ptms_list Cluster 5, Gene 6", {expect_equal(eu_ptms_list[[5]]$Gene.Name[6], "ACBD3 ubi K386")})

#sed_ptms_list
test_that("Testing sed_ptms_list Cluster 3, Gene 1", {expect_equal(sed_ptms_list[[3]]$Gene.Name[1], "AASDHPPT ack K151")})
test_that("Testing sed_ptms_list Cluster 4, Gene 1", {expect_equal(sed_ptms_list[[4]]$Gene.Name[1], "ABCA1 ubi K2023")})
test_that("Testing sed_ptms_list Cluster 5, Gene 9", {expect_equal(sed_ptms_list[[5]]$Gene.Name[9], "ABCC4 ubi K77")})

#sp_ptms_list
test_that("Testing sp_ptms_list Cluster 8, Gene 1", {expect_equal(sp_ptms_list[[8]]$Gene.Name[1], "ABCB6 ubi K482")})
test_that("Testing sp_ptms_list Cluster 13, Gene 1", {expect_equal(sp_ptms_list[[13]]$Gene.Name[1], "ABCC1 ubi K290")})
test_that("Testing sp_ptms_list Cluster 21, Gene 2", {expect_equal(sp_ptms_list[[21]]$Gene.Name[2], "ACKR3 ubi K362")})

#Clean Up
if(file.exists("noprint")) file.remove("noprint") #Clean up file created by sink()
