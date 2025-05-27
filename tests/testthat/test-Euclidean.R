#Set the seed (very important)
set.seed(1)

#Load Sample data
if(!exists("ptmtable")){
  load("~/CCCN_CFN_Tools/data/ptmtable.rda")
}

#Create sample data
sink("noprint") #Suppress print statements
testptmt <- EuclideanDistance(ptmtable)
sink()

#n is the tolerance (how close values of testptmt can be to their respective expected value)
n <- 1e-5

#Unit Testing for Euclidean Distance - IMPROVE ME
test_that("Testing EuclideanDistance Row 1, Col 1", {expect_equal(testptmt[1, 1], 17.41404, tolerance = n)})
test_that("Testing EuclideanDistance Row 2, Col 2", {expect_equal(testptmt[1, 2], -9.51064, tolerance = n)})
test_that("Testing EuclideanDistance Row 3, Col 3", {expect_equal(testptmt[1, 3], -1.58277, tolerance = n)})

if(file.exists("noprint")) file.remove("noprint") #Clean Up
