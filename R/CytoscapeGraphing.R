#' Graph Cluster Filtered Network
#'
#' Creates a cytoscape graph of the cluster filtered network. Ensure that you have the cytoscape app open and the RCy3 package downloaded and libraried.
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

  setNodeShapeDefault("VEE")                        # "RECTANGLE", "VEE", "OCTAGON", "ELLIPSE", "PARALLELOGRAM", "DIAMOND", "HEXAGON", "ROUND_RECTANGLE", "TRIANGLE"
  setNodeSizeDefault(50)                            # self explanatory
  setEdgeOpacityDefault(175)                        # 0 - 255
  setNodeLabelPositionDefault("C", "C", "c", 0, 0)  # What part of the node label is aligned "C", "NW", "N", "NE", "E", "SE", "S", "SW", "W"
                                                    # to what part of the node graphic       "C", "NW", "N", "NE", "E", "SE", "S", "SW", "W"
                                                    # "l", "r", "c"
                                                    # amount offset in the x direction
                                                    # amount offset in the y direction
  setBackgroundColorDefault('#EDC8FF')              # self explanatory
  setEdgeColorDefault('#899499')                    # self explanatory

  # mapVisualProperty(visual.property, table.column, mapping.type, table.column.values, visual.prop.values)
  # visual.properties can be found with ?mapVisualProperty (thank god) some useful ones are: "node label", "node fill color", "Edge Target Arrow Shape", many more
  # table.column shows which column of the table is being mapped (id for the node label) or analyzed (interaction for arrow shape) for this property
  # mapping.type is 'c' for continuous (you know what continuous means); 'd' for discrete (set determined vals as in arrow shape and node color); 'p' for passthrough (each is different? used for node label)
  # table.column.values determine which entries for that column give you which VVVV
  # visual.prop.values give the property type assigned to the specified ^^^^^

}
