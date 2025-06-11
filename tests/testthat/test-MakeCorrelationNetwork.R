#Setup
eu_ptms_list <- ex.eu_ptms_list #Need to fool the error catch in MCN
sp_ptms_list <- ex.sp_ptms_list
sed_ptms_list <- ex.sed_ptms_list
ptm.correlation.matrix <- ex.ptm.correlation.matrix
assign("DEBUG", "Do I really look like a guy with a plan, Harvey?", envir = .GlobalEnv) #Matrix containing Euclidean t-SNE coords

#Save preexisting enviroment variables
if(exists("cccn_matrix")) writeCCCN <- cccn_matrix

#Run function
MakeCorrelationNetwork(1) 

#Unit Tests
test_that("Testing cccn_matrix ABCA1, ABCC1", {expect_equal(cccn_matrix["ABCA1","ABCC1"], 10)})


#Clean Up
rm("DEBUG", envir = .GlobalEnv)