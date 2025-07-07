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
GraphCfn <- function(cfn, Network.title = "cfn", Network.collection = "cccn.cfn.tools",
                     background.color = '#fcf3cf', edge.label.color = '#17202a', edge.line.color = '#abb2b9', node.border.color = '#145a32', node.label.color = '#145a32', node.color = "#a9dfbf",
                     default.font = "Times New Roman", node.font.size = 12, edge.font.size = 8,
                     edge.line.style = 'SOLID', source.arrow = 'NONE', target.arrow = 'NONE', node.shape = "OCTAGON",
                     node.size = 50, edge.width = 2, border.width = 1,
                     edge.opacity = 175, edge.label.opacity = 255, border.opacity = 255, node.label.opacity = 255, node.fill.opacity = 255
                     ){

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

  # STANDARD: if want to change; do separately

  setEdgeSourceArrowColorDefault(edge.line.color)   # set color of arrow from source
  setEdgeTargetArrowColorDefault(edge.line.color)   # set color of arrow to target
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

  # CUSTOMIZED: user inputs can change the following:

  # colors
  setBackgroundColorDefault(background.color)       # set color of background
  setEdgeLabelColorDefault(edge.label.color)        # set color of edge label
  setEdgeColorDefault(edge.line.color)              # set color of edge
  setNodeBorderColorDefault(node.border.color)      # set color of node border
  setNodeLabelColorDefault(node.label.color)        # set color of node title
  setNodeColorDefault(node.color)                   # set interior color of node
  # fonts
  setEdgeFontFaceDefault(default.font)              # set font of edge
  setNodeFontFaceDefault(default.font)              # set font of node title (Initial Default UNKNOWN (not TNR or Arial))
  setEdgeFontSizeDefault(edge.font.size)            # set font size of edge (Initial Default 12)
  setNodeFontSizeDefault(node.font.size)            # set font size fo node title (Initial Default 12)
  # shape/style
  setEdgeLineStyleDefault(edge.line.style)          # "PARALLEL_LINES", "MARQUEE_EQUAL", "DOT", "EQUAL_DASH", "LONG_DASH", "CONTIGUOUS_ARROW", "MARQUEE_DASH", "DASH_DOT", "BACKWARD_SLASH", "FORWARD_SLASH", "VERTICAL_SLASH", "SOLID", "SEPARATE_ARROW", "MARQUEE_DASH_DOT", "ZIGZAG", "SINEWAVE"
  setEdgeSourceArrowShapeDefault(source.arrow)      # "DELTA", "DIAMOND", "OPEN_CIRCLE", "CIRCLE", "OPEN_HALF_CIRCLE", "CROSS_OPEN_DELTA", "DELTA_SHORT_1", "CROSS_DELTA", "OPEN_DELTA", "OPEN_DIAMOND"
  setEdgeTargetArrowShapeDefault(target.arrow)      # "DIAMOND_SHORT_1", "DELTA_SHORT_2", "OPEN_SQUARE", "NONE", "SQUARE", "DIAMOND_SHORT_2", "T", "HALF_BOTTOM", "HALF_TOP", "ARROW_SHORT", "HALF_CIRCLE"
  setNodeShapeDefault(node.shape)                   # "RECTANGLE", "VEE", "OCTAGON", "ELLIPSE", "PARALLELOGRAM", "DIAMOND", "HEXAGON", "ROUND_RECTANGLE", "TRIANGLE"
  # size
  setNodeSizeDefault(node.size)                     # set size of node; height and width assumed the same
  setEdgeLineWidthDefault(edge.width)               # set edge line thickness; Initial Default 2
  setNodeBorderWidthDefault(border.width)           # set border thickness; Initial Default 0
  # opacity
  setEdgeOpacityDefault(edge.opacity)               # set opacity of edge; 0 - 255 w 0 --> translucent
  setEdgeLabelOpacityDefault(edge.label.opacity)    # set opacity of edge label; 0 - 255 w 0 --> translucent
  setNodeBorderOpacityDefault(border.opacity)       # set opacity of border of node; 0 - 255 w 0 --> translucent
  setNodeFillOpacityDefault(node.fill.opacity)      # set opacity of interior color of node; 0 - 255 w 0 --> translucent
  setNodeLabelOpacityDefault(node.label.opacity)    # set opacity of title of node; 0 - 255 w 0 --> translucent

}
