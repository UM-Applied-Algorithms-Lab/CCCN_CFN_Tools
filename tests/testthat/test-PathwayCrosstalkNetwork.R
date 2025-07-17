#Load Sample Data
sink(file="noprint_123a") #This redirects all printing to a file named noprint_123a
PathwayCrosstalkNetwork(ex.common.clusters, ex.bioplanet, edgelist.name = "edgelist_123a")
sink()

#Unit tests for the structures in PCN
#Jaccard Weights
jaccardtest <- sapply(edgelist_123a[,3], as.double)
test_that("Row 1 Jaccard Weight", {expect_equal(jaccardtest[[1]], 0.00173310,tolerance=0.00001)})
test_that("Row 4 Jaccard Weight", {expect_equal(is.na(jaccardtest[[4]]), TRUE)})
test_that("Row 7 Jaccard Weight", {expect_equal(jaccardtest[[7]], 0.02597402,tolerance=0.00001)})

#CPE Evidence
cpetest <- sapply(edgelist_123a[,4], as.double)
test_that("Row 3 PTP interaction evidence", {expect_equal(cpetest[[3]], 1.851253,tolerance=0.00001)})
test_that("Row 4 PTP interaction evidence", {expect_equal(cpetest[[4]], 0.137931,tolerance=0.00001)})
test_that("Row 17 PTP interaction evidence", {expect_equal(is.na(cpetest[[17]]), TRUE)})

#Cleanup
rm("edgelist_123a", envir = .GlobalEnv) #remove the variable created
file.remove("noprint_123a") #remove the noprint file
file.remove("edgelist_123a.csv") #remove the edgelist file created
