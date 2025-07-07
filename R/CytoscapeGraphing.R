#' Graph Cluster Filtered Network
#'
#' Creates a cytoscape graph of the cluster filtered network. Ensure that you have the cytoscape app open and the RCy3 package downloaded and libraried.
#'
#' @param cfn A version of ppi.network with only the edges that exist in cccn.matrix and have non-zero weights
#' @param Network.title Desired title for the created Cytoscape Network
#' @param Network.collection Desired name for the collection created on Cytoscape in which the network will reside
#'
#' @return A cytoscape graph of the cluster filtered network
#' @export
#'
#' @examples
#' GraphCFN(ex.cfn)
GraphCfn <- function(cfn, Network.title = "cfn", Network.collection = "cccn.cfn.tools"){

  genes <- unique(c(cfn$Gene.1, cfn$Gene.2))

  cfn.edges <- data.frame(matrix(data = 0, nrow = length(rownames(cfn)), ncol = 4))
  cfn.nodes <- data.frame(matrix(data = 0, nrow = length(genes), ncol = 3))

  colnames(cfn.edges) <- c("source", "target", "interaction", "weight")
  colnames(cfn.nodes) <- c("id", "group", "score")

  cfn.edges$source <- cfn$Gene.1
  cfn.edges$target <- cfn$Gene.2
  cfn.edges$weight <- cfn$PPI.weight

  cfn.nodes$id <- genes

  createNetworkFromDataFrames(cfn.nodes, cfn.edges, title = Network.title, collection = Network.collection)

  # CUSTOMIZATION FROM HERE ON OUT

  setBackgroundColorDefault('#EDC8FF')              # self explanatory

  # edge stuff
  setEdgeFontFaceDefault()                          # set font of edge
  setEdgeFontSizeDefault()                          # set font size of edge
  setEdgeOpacityDefault(175)                        # 0 - 255
  setEdgeLabelColorDefault()                        #
  setEdgeLabelDefault()                             #
  setEdgeLabelOpacityDefault()                      #
  setEdgeLineStyleDefault()                         #
  setEdgeLineWidthDefault()                         #
  setEdgeOpacityDefault()                           #
  setEdgeColorDefault('#899499')                    # self explanatory
  # arrow stuff
  setEdgeSourceArrowColorDefault()                  #
  setEdgeSourceArrowShapeDefault()                  #
  setEdgeTargetArrowColorDefault()                  #
  setEdgeTargetArrowShapeDefault()                  #

  # node stuff
  setNodeShapeDefault("VEE")                        # "RECTANGLE", "VEE", "OCTAGON", "ELLIPSE", "PARALLELOGRAM", "DIAMOND", "HEXAGON", "ROUND_RECTANGLE", "TRIANGLE"
  setNodeSizeDefault(50)                            # self explanatory
  # border stuff
  setNodeBorderColorDefault()                       #
  setNodeBorderOpacityDefault()                     #
  setNodeBorderWidthDefault()                       #
  # else
  setNodeColorDefault()                             #
  setNodeFillOpacityDefault()                       #
  setNodeFontFaceDefault()                          #
  setNodeFontSizeDefault()                          #
  setNodeHeightDefault()                            #
  setNodeLabelColorDefault()                        #
  setNodeLabelOpacityDefault()                      #
  setNodeWidthDefault()                             #
  setNodeLabelPositionDefault("C", "C", "c", 0, 0)  # What part of the node label is aligned "C", "NW", "N", "NE", "E", "SE", "S", "SW", "W"
                                                    # to what part of the node graphic       "C", "NW", "N", "NE", "E", "SE", "S", "SW", "W"
                                                    # "l", "r", "c"
                                                    # amount offset in the x direction
                                                    # amount offset in the y direction

  # mapVisualProperty(visual.property, table.column, mapping.type, table.column.values, visual.prop.values)
  # visual.properties can be found with ?mapVisualProperty (thank god) some useful ones are: "node label", "node fill color", "Edge Target Arrow Shape", many more
  # table.column shows which column of the table is being mapped (id for the node label) or analyzed (interaction for arrow shape) for this property
  # mapping.type is 'c' for continuous (you know what continuous means); 'd' for discrete (set determined vals as in arrow shape and node color); 'p' for passthrough (each is different? used for node label)
  # table.column.values determine which entries for that column give you which VVVV
  # visual.prop.values give the property type assigned to the specified ^^^^^

}
