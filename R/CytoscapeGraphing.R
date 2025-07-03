#' Graph Cluster Filtered Network
#'
#' Creates a cytoscape graph of the cluster filtered network.
#'
#' @param cfn A version of ppi.network with only the edges that exist in cccn.matrix and have non-zero weights
#'
#' @return A cytoscape graph of the cluster filtered network
#' @export
#'
#' @examples
#' gmpath <- "genemania-interactions.txt"
#' pack <- "cccn.cfn.tools"
#' gmfile <- system.file("genemania", gmpath, package = pack, mustWork = TRUE)
#' cccn.cfn.tools:::ex.FindPPIEdges(ex.cccn.matrix, ppi.network.name = "ex.ppi.network")
GraphCfn <- function(cfn){

}
