#Load Sample Data
sink(file="noprint_123a") #This redirects all printing to a file named noprint_123a
PathwayCrosstalkNetwork(ex.bioplanet, ex.common.clusters, edgelist.name = "edgelist_123a")
sink()

#Unit tests for the structures in PCN
#Jaccard Weights
test_that("Row 1 Jaccard Weight", {expect_equal(edgelist_123a[[1,3]], 0.00173310,tolerance=0.00001)})
test_that("Row 4 Jaccard Weight", {expect_equal(edgelist_123a[[4,3]], NA,tolerance=0.00001)})
test_that("Row 7 Jaccard Weight", {expect_equal(edgelist_123a[[7,3]], 0.02597402,tolerance=0.00001)})

#CPE Evidence
test_that("Row 3 PTP interaction evidence", {expect_equal(edgelist_123a[[3,4]], 0.025588660,tolerance=0.00001)})
test_that("Row 4 PTP interaction evidence", {expect_equal(edgelist_123a[[4,4]], 0.024672242,tolerance=0.00001)})
test_that("Row 17 PTP interaction evidence", {expect_equal(edgelist_123a[[17,4]], NA,tolerance=0.00001)})

#Cleanup
rm("edgelist_123a", envir = .GlobalEnv) #remove the variable created
file.remove("noprint_123a") #remove the noprint file
file.remove("edgelist_123a.csv") #remove the edgelist file created
