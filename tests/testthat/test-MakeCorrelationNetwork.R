#Run function
MakeCorrelationNetwork(ex.clusters.list, ex.ptm.correlation.matrix, 1, clusters.name = "common_123a", cccn.name = "cccn_123a") 

#Test common clusters
test_that("Cluster 1, PTM 2", {expect_equal(common.clusters$`Cluster 1`[[2]], "ABCC1 ubi K1502")})
test_that("Cluster 6, PTM 1", {expect_equal(common.clusters$`Cluster 6`[[1]], "ABCA3 p S683")})
test_that("Cluster 7, PTM 1", {expect_equal(common.clusters$`Cluster 7`[[1]], "ABL1 p Y393; ABL2 p Y439")})
test_that("Cluster 8, PTM 2", {expect_equal(common.clusters$`Cluster 8`[[2]], "ABLIM1 p Y373")})

#Test cccn_matrix for correctness
test_that("Testing cccn_matrix ABCA1, ABCC1", {expect_equal(cccn.matrix["ABCA1","ABCC1"], 10)})
test_that("Testing cccn_matrix ACKR3, ABCC1", {expect_equal(cccn.matrix["ACKR3","ABCC1"], 9.6)})
test_that("Testing cccn_matrix ABLIM1, ABHD11", {expect_equal(cccn.matrix["ABLIM1","ABHD11"], 5.8)})
test_that("Testing cccn_matrix AATF, ABCD3", {expect_equal(cccn.matrix["AATF","ABCD3"], 1)})

#Clean Up
rm("common_123a", envir = .GlobalEnv)
rm("cccn_123a", envir = .GlobalEnv)