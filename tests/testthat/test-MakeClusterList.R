#Setup
set.seed(1) #Set the seed (very important)

#Load Sample data
sink("noprint")              #Suppress print statements from function
temp(ex.ptmtable) #Create sample data - #BUG - writes 'species scores not available' (dont worry about this for now)
sink()

#Unit Tests for the global variables (lists of dataframes) created by MakeClusterList. Check PTMs Names.
#eu_ptms_list
test_that("Testing eu_ptms_list Cluster 1, PTM 2", {expect_equal(eu_ptms_list[[1]]$PTM.Name[2], "ABCA1 ubi K2023")})
test_that("Testing eu_ptms_list Cluster 3, PTM 4", {expect_equal(eu_ptms_list[[3]]$PTM.Name[4], "ABCB6 ubi K482")})
test_that("Testing eu_ptms_list Cluster 5, PTM 6", {expect_equal(eu_ptms_list[[5]]$PTM.Name[6], "ACBD3 ubi K386")})

#sed_ptms_list
test_that("Testing sed_ptms_list Cluster 3, PTM 1", {expect_equal(sed_ptms_list[[3]]$PTM.Name[1], "AASDHPPT ack K151")})
test_that("Testing sed_ptms_list Cluster 4, PTM 1", {expect_equal(sed_ptms_list[[4]]$PTM.Name[1], "ABCB6 ubi K726")})

#sp_ptms_list
test_that("Testing sp_ptms_list Cluster 8, PTM 1", {expect_equal(sp_ptms_list[[8]]$PTM.Name[1], "ABCB6 ubi K482")})
test_that("Testing sp_ptms_list Cluster 13, PTM 1", {expect_equal(sp_ptms_list[[13]]$PTM.Name[1], "ABCC1 ubi K290")})
test_that("Testing sp_ptms_list Cluster 21, PTM 2", {expect_equal(sp_ptms_list[[21]]$PTM.Name[2], "ACKR3 ubi K362")})

#Cleanup
if(file.exists("noprint")) file.remove("noprint") #Clean up file created by sink()