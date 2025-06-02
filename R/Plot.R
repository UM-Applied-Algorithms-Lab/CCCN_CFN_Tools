#' A helper function that plots data from the outputs of various different functions
#'
#' @param fun If you pass in the name of a specific function you want to plot, then it will skip the asking step. A dot is default and will run the asking step
#'
#' @examples
#' Plot("MakeClusterList")
Plot <- function(fun = "."){

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

  if(fun=="."){
    #Find what has been run


  }

}
