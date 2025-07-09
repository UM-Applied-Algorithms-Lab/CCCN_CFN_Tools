#' Graph Cluster Filtered Network
#'
#' Creates a cytoscape graph of the cluster filtered network. Ensure that you have the cytoscape app open and the RCy3 package downloaded and libraried.
#'
#' If error occurs, run cytoscapePing() to ensure connection to the Cytoscape Interface and try again.
#'
#' @param cfn A version of ppi.network with only the edges that exist in cccn.matrix and have non-zero weights
#' @param Network.title Desired title for the created Cytoscape Network
#' @param Network.collection Desired name for the collection created on Cytoscape in which the network will reside
#'
#' @param background.color Hex code of background color of graph; defaults to '#fcf3cf'
#' @param edge.label.color Hex code of edge label color of graph; defaults to '#17202a'
#' @param edge.line.color Hex code of edge line color of graph; defaults to '#abb2b9'
#' @param node.border.color Hex code of node border color of graph; defaults to '#145a32'
#' @param node.label.color Hex code of node label color of graph; defaults to '#145a32'
#' @param node.fill.color Hex code of node fill color of graph; defaults to '#a9dfbf'
#'
#' @param default.font Font style of edge and node names; defaults to "Times New Roman"
#' @param node.font.size Font size of the node name; defaults to 12
#' @param edge.font.size Font size of the edge name; defaults to 8
#'
#' @param edge.line.style Type of edge style; defaults to "SOLID"; options include: "PARALLEL_LINES", "MARQUEE_EQUAL", "DOT", "EQUAL_DASH", "LONG_DASH", "CONTIGUOUS_ARROW", "MARQUEE_DASH", "DASH_DOT", "BACKWARD_SLASH", "FORWARD_SLASH", "VERTICAL_SLASH", "SOLID", "SEPARATE_ARROW", "MARQUEE_DASH_DOT", "ZIGZAG", "SINEWAVE"
#' @param source.arrow Type of arrow coming from the source gene; defaults to "NONE"; options include: "DELTA", "DIAMOND", "OPEN_CIRCLE", "CIRCLE", "OPEN_HALF_CIRCLE", "CROSS_OPEN_DELTA", "DELTA_SHORT_1", "CROSS_DELTA", "OPEN_DELTA", "OPEN_DIAMOND", "DIAMOND_SHORT_1", "DELTA_SHORT_2", "OPEN_SQUARE", "NONE", "SQUARE", "DIAMOND_SHORT_2", "T", "HALF_BOTTOM", "HALF_TOP", "ARROW_SHORT", "HALF_CIRCLE"
#' @param target.arrow Type of arrow going to the target gene; defaults to "NONE"; options include: "DELTA", "DIAMOND", "OPEN_CIRCLE", "CIRCLE", "OPEN_HALF_CIRCLE", "CROSS_OPEN_DELTA", "DELTA_SHORT_1", "CROSS_DELTA", "OPEN_DELTA", "OPEN_DIAMOND", "DIAMOND_SHORT_1", "DELTA_SHORT_2", "OPEN_SQUARE", "NONE", "SQUARE", "DIAMOND_SHORT_2", "T", "HALF_BOTTOM", "HALF_TOP", "ARROW_SHORT", "HALF_CIRCLE"
#' @param node.shape Shape of the node; defaults to "OCTAGON"; options include: "RECTANGLE", "VEE", "OCTAGON", "ELLIPSE", "PARALLELOGRAM", "DIAMOND", "HEXAGON", "ROUND_RECTANGLE", "TRIANGLE"
#'
#' @param node.size Size of the node; defaults to 50. PLEASE NOTE: width and height can be changed independently using RCy3 directly (first run lockNodeDimensions(FALSE) and then setNodeWidthDefault() and setNodeHeightDefault())
#' @param edge.width Width of the edge line; defaults to 2
#' @param border.width Width of the border of the node; defaults to 1
#'
#' @param edge.opacity Opacity of the edge line on a scale of 0 - 255 with 0 being transparent; defaults to 175
#' @param edge.label.opacity Opacity of the edge label on a scale of 0 - 255 with 0 being transparent; defaults to 255
#' @param border.opacity Opacity of the node border on a scale of 0 - 255 with 0 being transparent; defaults to 255
#' @param node.label.opacity Opacity of the node label on a scale of 0 - 255 with 0 being transparent; defaults to 255
#' @param node.fill.opacity Opacity of the node fill on a scale of 0 - 255 with 0 being transparent; defaults to 255
#'
#' @return A cytoscape graph of the cluster filtered network
#' @export
#'
#' @examples
#' # GraphCFN(ex.cfn)
#' # See vignette for default graph
GraphCfn <- function(cfn, Network.title = "cfn", Network.collection = "cccn.cfn.tools",
                     background.color = '#fcf3cf', edge.label.color = '#17202a', edge.line.color = '#abb2b9', node.border.color = '#145a32', node.label.color = '#145a32', node.fill.color = '#a9dfbf',
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
    cytoscapePing()                           # ensure connection?
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

  # mapVisualProperty("Edge Target Arrow Shape", "interaction", 'd',
  #                   c("Physical Interactions", "Pathways", "database", "database_transferred", "experimental", "experimental_transferred"),
  #                   c("SQUARE", "HALF_CIRCLE", "CROSS_DELTA", "CROSS_DELTA", "ARROW_SHORT", "ARROW_SHORT"))

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
  setNodeLabelColorDefault(node.label.color)        # set color of node name
  setNodeColorDefault(node.fill.color)              # set interior color of node
  # fonts
  setEdgeFontFaceDefault(default.font)              # set font of edge
  setNodeFontFaceDefault(default.font)              # set font of node name (Initial Default UNKNOWN (not TNR or Arial))
  setEdgeFontSizeDefault(edge.font.size)            # set font size of edge (Initial Default 12)
  setNodeFontSizeDefault(node.font.size)            # set font size of node name (Initial Default 12)
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
  setNodeLabelOpacityDefault(node.label.opacity)    # set opacity of name of node; 0 - 255 w 0 --> translucent

}
