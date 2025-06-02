#This file might be completely obsolete but will probably copy and paste functions from it?

#' A helper function that plots data from the outputs of various different functions
#'
#' @param fun If you pass in the name of a specific function you want to plot, then it will skip the asking step. A dot is default and will run the asking step
#'
#' @examples
#' Plot("MakeClusterList")
Plot <- function(dir = "plots", fun = "."){

  #Plot MakeClusterList
  PlotClusters <- function(){
    print("ran!") #TO DO
  }

  #Plot MakeCorrelationNetwork
  PlotCorrelationNetwork <- function(){
    print("run!") #TO DO
  }

  #Plot FindPPIEdges
  PlotPPI <- function(){
    print("runned!") #TO DO
  }

  #Plot ClusterFilteredNetwork
  PlotClusterFilteredNetwork <- function(){
    print("ranned!") #TO DO
  }

  #Init
  ranMCL <- FALSE
  ranClusterFN <- FALSE
  ranPPI <- FALSE
  ranMCN <- FALSE

  #Find what has been run: (Please for the love of good I think I'd be so much easier)
  if((exists("ptm.correlation.matrix"))) ranMCL = TRUE
  if((exists("Network"))) ranMCN                = TRUE #Does NOT work since MCN does not assign anything to global namespace
  if((exists("combined_ppi_network"))) ranPPI   = TRUE
  if((exists("???"))) ranClusterFN              = TRUE

  #Plot help?
  if(fun == "."){
    if(!ranMCL || !ranMCN || !ranPPI || !ranClusterFN) stop("Please run a function before running plot. See vignette for a list of functions")

    #If at least 1 function has been ran

    if(ranMCL) print("")
  }



}
