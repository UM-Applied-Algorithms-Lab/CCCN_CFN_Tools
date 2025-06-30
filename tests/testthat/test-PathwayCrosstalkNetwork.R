#Load Sample Data

PathwayCrosstalkNetwork(ex.bioplanet, ex.common.clusters, edgelist.name = "edgelist_123a")

#Unit tests for the structures in PCN
#Jaccard Weights
writeLines(""); writeLines("Testing Edgelist Jaccard weights")
test_that("Row 10 Jaccard Weight", {expect_equal(edgelist_123a$Jaccard_weight[[10]], 0.1428571,tolerance=0.00001)})
test_that("Row 20 Jaccard Weight", {expect_equal(edgelist_123a$Jaccard_weight[[20]], 0,tolerance=0.00001)})
test_that("Row 66 Jaccard Weight", {expect_equal(edgelist_123a$Jaccard_weight[[66]], 0.4285714,tolerance=0.00001)})

#CPE Evidence
test_that("Row 15 PTP interaction evidence", {expect_equal(edgelist_123a$PTP_inter_evidence[[15]], 0.1666667,tolerance=0.00001)})
test_that("Row 35 PTP interaction evidence", {expect_equal(edgelist_123a$PTP_inter_evidence[[35]], 0.2083333,tolerance=0.00001)})
test_that("Row 56 PTP interaction evidence", {expect_equal(edgelist_123a$PTP_inter_evidence[[56]], 1,tolerance=0.00001)})

#Cleanup
rm("edgelist_123a", envir = .GlobalEnv) #remove the variable created
file.remove("edgelist_123a.csv") #remove the edgelist file created