#Setup
set.seed(1) #Set the seed (very important)

#Load Sample data
MakeClusterList(ex.ptmtable, name.columns = 1:3, correlation.matrix.name = "matrix_123a", clusters.list.name = "clusters_123a", tsne.coords.name = "tsne_123a", toolong = 3.5) #Create sample data - #BUG - writes 'species scores not available' (dont worry about this for now)

#Unit Tests for the global variables (lists of dataframes) created by MakeClusterList. Check PTMs Names.
#eu_ptms_list
test_that("Cluster 2, PTM 2", {expect_equal(clusters_123a$Euclidean$`2`$PTMnames[2], "ABLIM1 p Y410")})
test_that("Cluster 3, PTM 4", {expect_equal(clusters_123a$Euclidean$`3`$PTMnames[4], "HNRNPU ack K352")})
test_that("Cluster 5, PTM 6", {expect_equal(clusters_123a$Euclidean$`5`$PTMnames[6], "HSPA1B ubi K597")})

#sed_ptms_list
test_that("Testing Cluster 3, PTM 1", {expect_equal(clusters_123a$SED$`3`$PTMnames[1], "ABLIM1 p Y383")})
test_that("Testing Cluster 4, PTM 1", {expect_equal(clusters_123a$SED$`4`$PTMnames[1], "ABLIM1 p Y406")})

#sp_ptms_list
test_that("Testing Cluster 8, PTM 1", {expect_equal(clusters_123a$Spearman$`8`$PTMnames[1], "ABLIM1 ubi K8")})
test_that("Testing Cluster 13, PTM 1", {expect_equal(clusters_123a$Spearman$`13`$PTMnames[1], "APH1A ubi K95")})
test_that("Testing Cluster 21, PTM 2", {expect_equal(clusters_123a$Spearman$`21`$PTMnames[2], "ATP6V0D1 p Y347")})

#Cleanup
rm(matrix_123a, envir = .GlobalEnv)
rm(clusters_123a, envir = .GlobalEnv)
rm(tsne_123a, envir = .GlobalEnv)
