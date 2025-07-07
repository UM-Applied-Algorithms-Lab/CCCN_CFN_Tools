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

  if(!(exists("RCy3"))){                      # check if RCy3 is libraried
    if(system.file(package="RCy3") == ""){    # check if RCy3 is installed at all
      BiocManager::install(version='devel')   # The following initializes usage of Bioc devel
      BiocManager::install("RCy3")            # install!
    }
    library("RCy3")                           # library
  }





  # ACTUAL CODE AND DATA PROCESSING

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
  setEdgeLabelColorDefault()                        # self explanatory
  setEdgeLabelOpacityDefault()                      # 0 - 255
  setEdgeLineStyleDefault()                         # "PARALLEL_LINES", "MARQUEE_EQUAL", "DOT", "EQUAL_DASH", "LONG_DASH", "CONTIGUOUS_ARROW", "MARQUEE_DASH", "DASH_DOT", "BACKWARD_SLASH", "FORWARD_SLASH", "VERTICAL_SLASH", "SOLID", "SEPARATE_ARROW", "MARQUEE_DASH_DOT", "ZIGZAG", "SINEWAVE"
  setEdgeLineWidthDefault()                         # Initial Default 2
  setEdgeOpacityDefault()                           # 0 - 255
  setEdgeColorDefault('#899499')                    # self explanatory
  # arrow stuff
  setEdgeSourceArrowColorDefault('#899499')         # self explanatory
  setEdgeTargetArrowColorDefault('#899499')         # self explanatory
  setEdgeSourceArrowShapeDefault()                  # "DELTA", "DIAMOND", "OPEN_CIRCLE", "CIRCLE", "OPEN_HALF_CIRCLE", "CROSS_OPEN_DELTA", "DELTA_SHORT_1", "CROSS_DELTA", "OPEN_DELTA", "OPEN_DIAMOND"
  setEdgeTargetArrowShapeDefault()                  # "DIAMOND_SHORT_1", "DELTA_SHORT_2", "OPEN_SQUARE", "NONE", "SQUARE", "DIAMOND_SHORT_2", "T", "HALF_BOTTOM", "HALF_TOP", "ARROW_SHORT", "HALF_CIRCLE"

  # node stuff
  setNodeShapeDefault("VEE")                        # "RECTANGLE", "VEE", "OCTAGON", "ELLIPSE", "PARALLELOGRAM", "DIAMOND", "HEXAGON", "ROUND_RECTANGLE", "TRIANGLE"
  setNodeSizeDefault(50)                            # self explanatory
  # border stuff
  setNodeBorderColorDefault('#899499')              # self explanatory
  setNodeBorderOpacityDefault(255)                  # 0 - 255
  setNodeBorderWidthDefault(0)                      # Initial Default 0
  # else
  setNodeColorDefault("#89D0F5")                    # self explanatory
  setNodeFillOpacityDefault()                       # 0 - 255
  setNodeFontFaceDefault()                          #
  setNodeFontSizeDefault()                          #
  setNodeSizeDefault()                              # Initial Default 50
  setNodeLabelColorDefault('#899499')               # self explanatory
  setNodeLabelOpacityDefault()                      # 0 - 255
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
