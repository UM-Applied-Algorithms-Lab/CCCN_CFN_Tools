#Run function
MakeCorrelationNetwork(ex.common.clusters, ex.ptm.correlation.matrix, ptm.cccn.name = "ptms_123a", gene.cccn.name = "gene_123a") 


#Test ptms, since PTMs are just a filtered version of the correlation matrix and thats tested already, only 2 tests
test_that("Testing ptms matrix HSPA1B ubi K597, PKP4 p Y1156", {expect_equal(ptms_123a["HSPA1B ubi K597","PKP4 p Y1156"], 0.4285714, tolerance=0.00001)})


#Test gene matrix for correctness
test_that("Testing gene network diagionals are 0", {expect_equal(gene_123a[43, 43], 0)})
test_that("Testing gene network CTTN, SEC16A", {expect_equal(gene_123a["CTTN", "SEC16A"], 60.94232, tolerance=0.00001)})
test_that("Testing gene network ARAP2, OAT", {expect_equal(gene_123a["ARAP2", "OAT"], 0)})
test_that("Testing gene network FKBP4, NSFL1C", {expect_equal(gene_123a["FKBP4", "NSFL1C"], 1.628571, tolerance=0.00001)})


#Clean Up
rm("ptms_123a", envir = .GlobalEnv)
rm("gene_123a", envir = .GlobalEnv)