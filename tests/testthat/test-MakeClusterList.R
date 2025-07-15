#Setup
set.seed(1) #Set the seed (very important)

#Load Sample data
MakeClusterList(ex.ptmtable, name.columns = 1:3, "matrix_123a", "clusters_123a", "tsne_123a") #Create sample data - #BUG - writes 'species scores not available' (dont worry about this for now)

#Unit Tests for the global variables (lists of dataframes) created by MakeClusterList. Check PTMs Names.
#eu_ptms_list
test_that("Cluster 2, PTM 2", {expect_equal(clusters_123a$Euclidean$`2`$PTM.Name[2], "AATF ack K145")})
test_that("Cluster 3, PTM 4", {expect_equal(clusters_123a$Euclidean$`3`$PTM.Name[4], "ABCB6 ubi K482")})
test_that("Cluster 5, PTM 6", {expect_equal(clusters_123a$Euclidean$`5`$PTM.Name[6], "ACBD3 ubi K386")})

#sed_ptms_list
test_that("Testing Cluster 3, PTM 1", {expect_equal(clusters_123a$SED$`3`$PTM.Name[1], "AASDHPPT ack K151")})
test_that("Testing Cluster 4, PTM 1", {expect_equal(clusters_123a$SED$`4`$PTM.Name[1], "ABCB6 ubi K726")})

#sp_ptms_list
test_that("Testing Cluster 8, PTM 1", {expect_equal(clusters_123a$Spearman$`8`$PTM.Name[1], "ABCB6 ubi K482")})
test_that("Testing Cluster 13, PTM 1", {expect_equal(clusters_123a$Spearman$`13`$PTM.Name[1], "ABCC1 ubi K290")})
test_that("Testing Cluster 21, PTM 2", {expect_equal(clusters_123a$Spearman$`21`$PTM.Name[2], "ACKR3 ubi K362")})

#Cleanup
rm(matrix_123a, envir = .GlobalEnv)
rm(clusters_123a, envir = .GlobalEnv)
rm(tsne_123a, envir = .GlobalEnv)
