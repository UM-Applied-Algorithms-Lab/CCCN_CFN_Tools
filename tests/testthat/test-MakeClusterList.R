#Setup
set.seed(1) #Set the seed (very important)

#Load Sample data
MakeClusterList(ex.ptmtable, correlation.matrix.name = "matrix_123a", clusters.list.name = "clusters_123a", tsne.coords.name = "tsne_123a", common.clusters.name = "common_123a", keeplength = 2, toolong = 3.5) #Create sample data - #BUG - writes 'species scores not available' (dont worry about this for now)



#Unit tests for the values of the correlation matrix
test_that("Diagionals are NA", {expect_equal(is.na(matrix_123a[10,10]), TRUE)})
test_that("APP p Y757 correlation with APP ubi K763", {expect_equal(matrix_123a[18, 20], 0.02857143,tolerance=0.00001)})
test_that("CTTN p Y334 correlation with FRS2 p Y306", {expect_equal(matrix_123a[64, 128], 0.912892,tolerance=0.00001)})



#Unit Tests for the global variables (lists of dataframes) created by MakeClusterList. Check PTMs Names.
#eu_ptms_list
test_that("Cluster 2, PTM 2", {expect_equal(clusters_123a$Euclidean$`2`$PTMnames[2], "ABLIM1 p Y410")})
test_that("Cluster 3, PTM 4", {expect_equal(clusters_123a$Euclidean$`3`$PTMnames[4], "HNRNPU ack K352")})
test_that("Cluster 5, PTM 6", {expect_equal(clusters_123a$Euclidean$`5`$PTMnames[6], "EZR ack K60")})

#sed_ptms_list
test_that("Cluster 3, PTM 1", {expect_equal(clusters_123a$SED$`3`$PTMnames[1], "ABLIM1 p Y406")})
test_that("Cluster 4, PTM 1", {expect_equal(clusters_123a$SED$`4`$PTMnames[1], "ABLIM1 p Y461")})

#sp_ptms_list
test_that("Cluster 8, PTM 1", {expect_equal(clusters_123a$Spearman$`8`$PTMnames[1], "APH1A p Y256")})
test_that("Cluster 13, PTM 1", {expect_equal(clusters_123a$Spearman$`13`$PTMnames[1], "APP ubi K763")})
test_that("Cluster 21, PTM 2", {expect_equal(clusters_123a$Spearman$`21`$PTMnames[2], "CTTN p Y199")})



#Unit tests for common clusters
test_that("SLC12A4 ubi K979 in cluster 4", {expect_equal("SLC12A4 ubi K979" %in% common_123a[[4]], TRUE)})
test_that("PKP4 p Y369 in cluster 10", {expect_equal("PKP4 p Y369" %in% common_123a[[10]], TRUE)})
test_that("OAT p Y126 in cluster 6", {expect_equal("OAT p Y126" %in% common_123a[[6]], TRUE)})




#Cleanup
rm(matrix_123a, envir = .GlobalEnv)
rm(clusters_123a, envir = .GlobalEnv)
rm(tsne_123a, envir = .GlobalEnv)
