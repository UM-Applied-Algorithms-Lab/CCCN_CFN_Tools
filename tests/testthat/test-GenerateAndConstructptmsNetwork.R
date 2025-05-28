#Setup
set.seed(1)                       #Set the seed (very important)

#Load Sample data
load("~/CCCN_CFN_Tools/data/ptmtable.rda")
par(mar=c(1,1,1,1)) #Resize plot

#Run function
sink("noprint")
GenerateAndConstructptmsNetwork(ptmtable, 0)
sink()

#UNIT TESTS
#Are plots being made correctly?
test_that("Check to see if the directory is made", {expect_equal(dir.exists("plots"), TRUE)})

#Check output Data
test_that("Check output data", {expect_equal(ptmtable$H3122SEPTM.C1[9], 20.29094)})
test_that("Check output data", {expect_equal(ptmtable$H3122SEPTM.C1[10], 18.22727)})
test_that("Check output data", {expect_equal(ptmtable$H2286_DMSO[93], 25.16108)})




#Cleanup
if(file.exists("noprint")) file.remove("noprint") #Clean up file created by sink()
unlink("plots", recursive=TRUE)
