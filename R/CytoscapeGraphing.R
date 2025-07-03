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
#' GraphCFN(ex.cfn)
GraphCfn <- function(cfn, network.name = "cfn.cytoscape"){

  genes <- unique(c(cfn$Gene.1, cfn$Gene.2))

  cfn.edges <- data.frame(matrix(data = 0, nrow = length(rownames(cfn)), ncol = 4))
  cfn.nodes <- data.frame(matrix(data = 0, nrow = length(genes), ncol = 3))

  colnames(cfn.edges) <- c("source", "target", "interaction", "weight")
  colnames(cfn.nodes) <- c("id", "group", "score")

  cfn.edges$source <- cfn$Gene.1
  cfn.edges$target <- cfn$Gene.2
  cfn.edges$weight <- cfn$PPI.weight

  cfn.nodes$id <- genes

  createNetworkFromDataFrames(cfn.nodes, cfn.edges, title = "cfn-test-network", collection = "cfn-testing")

}
