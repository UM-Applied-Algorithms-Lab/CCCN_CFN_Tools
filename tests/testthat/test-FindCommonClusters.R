#Setup
if(!exists("Sample3_ptm_list.rda")){
  load("~/CCCN_CFN_Tools/tests/testthat/Sample1_ptm_list.rda")
  load("~/CCCN_CFN_Tools/tests/testthat/Sample2_ptm_list.rda")
  load("~/CCCN_CFN_Tools/tests/testthat/Sample3_ptm_list.rda")
}

#Sample data to test with
common <- FindCommonClusters(Sample1_ptm_list, Sample2_ptm_list, Sample3_ptm_list, 2)

#UNIT TESTS FOR 3 FULL LISTS
test_that("Testing common cluster 1", {expect_equal(common[[1]], c("Apple", "Avocado", "Banana", "Blackberry"))})
test_that("Testing common cluster 2", {expect_equal(common[[2]], c("Cherry", "Citron", "Date"))})
test_that("Testing common cluster 5", {expect_equal(common[[5]], c("Starfruit", "Strawberry", "Tayberry"))})

#UNIT TEST FOR 2 FULL LISTS
common2 <- FindCommonClusters(Sample1_ptm_list, Sample2_ptm_list, list(), 2)
test_that("Testing common cluster 1", {expect_equal(common[[1]], c("Apple", "Avocado", "Banana", "Blackberry"))})


