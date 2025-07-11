#' Graph Cluster Filtered Network
#'
#' Creates a cytoscape graph of the cluster filtered network. Ensure that you have the cytoscape app open and the RCy3 package downloaded and libraried.
#'
#' If error occurs, run cytoscapePing() to ensure connection to the Cytoscape Interface and try again.
#'
#' @param cfn A version of ppi.network with only the edges that exist in cccn.matrix and have non-zero weights
#' @param ptmtable A dataset for post-translational modifications. Formatted with numbered rows, and the first column containing PTM names. The rest of the column names should be drugs. Values are numeric values that represent how much the PTM has reacted to the drug.
#' @param funckey A table graphing gene names to type of protein; defaults to internal database at cccn.cfn.tools::ex.funckey
#' @param Network.title Desired title for the created Cytoscape Network; defaults to "cfn"
#' @param Network.collection Desired name for the collection created on Cytoscape in which the network will reside; defaults to "cccn.cfn.tools"
#' @param visual.style.name Desired name for the visual style created on Cytoscape; defaults to "cccn.cfn.tools.style"
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
GraphCfn <- function(cfn, ptmtable, funckey = cccn.cfn.tools::ex.funckey, Network.title = "cfn", Network.collection = "cccn.cfn.tools", visual.style.name = "cccn.cfn.tools.style",
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

  ptmnew <- ptmtable[, c("PTM", "HCC4006_Erlotinib")]
  ptmnew$PTM <- sapply(ptmtable$PTM, function(x) strsplit(x, split = " ")[[1]][1])
  colnames(ptmnew) <- c("PTM", "score")
  ptmnew$score <- sapply(1:length(rownames(ptmtable)), function(i) min(as.numeric(ptmtable[-1][i, ]), na.rm = TRUE))

  cfn.edges <- data.frame(matrix(data = 0, nrow = length(rownames(cfn)), ncol = 4), stringsAsFactors = FALSE)
  cfn.nodes <- data.frame(matrix(data = 0, nrow = length(genes), ncol = 3), stringsAsFactors = FALSE)

  colnames(cfn.edges) <- c("source", "target", "interaction", "weight")
  colnames(cfn.nodes) <- c("id", "node.type", "score")

  cfn.edges$source <- cfn$Gene.1
  cfn.edges$target <- cfn$Gene.2
  cfn.edges$interaction <- cfn$Interaction
  cfn.edges$weight <- cfn$PPI.weight

  cfn.nodes$id <- genes
  cfn.nodes$node.type <- as.character(sapply(cfn.nodes$id, function(x) funckey$nodeType[which(funckey$Gene.Name == x)]))
  cfn.nodes$score <- sapply(cfn.nodes$id, function(x) sum(ptmnew$score[which(ptmnew$PTM == x)]))

  createNetworkFromDataFrames(cfn.nodes, cfn.edges, title = Network.title, collection = Network.collection)

  copyVisualStyle("default", visual.style.name)



  # CUSTOMIZATION FROM HERE ON OUT

  # STANDARD: if want to change; do separately

  setEdgeSourceArrowColorDefault(edge.line.color, visual.style.name)   # set color of arrow from source
  setEdgeTargetArrowColorDefault(edge.line.color, visual.style.name)   # set color of arrow to target
  setNodeLabelPositionDefault("C", "C", "c", 0, 0, visual.style.name)  # What part of the node label is aligned "C", "NW", "N", "NE", "E", "SE", "S", "SW", "W"
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
  # mapping.type is 'c' for continuous (you know what continuous means); 'd' for discrete (set determined vals as in arrow shape and node color); 'p' for passthrough (represented by the val in the table)
  # table.column.values determine which entries for that column give you which VVVV
  # visual.prop.values give the property type assigned to the specified ^^^^^

  # CUSTOMIZED: user inputs can change the following:

  # colors
  setBackgroundColorDefault(background.color, visual.style.name)       # set color of background
  setEdgeLabelColorDefault(edge.label.color, visual.style.name)        # set color of edge label
  setEdgeColorDefault(edge.line.color, visual.style.name)              # set color of edge
  setNodeBorderColorDefault(node.border.color, visual.style.name)      # set color of node border
  setNodeLabelColorDefault(node.label.color, visual.style.name)        # set color of node name
  setNodeColorDefault(node.fill.color, visual.style.name)              # set interior color of node
  # fonts
  setEdgeFontFaceDefault(default.font, visual.style.name)              # set font of edge
  setNodeFontFaceDefault(default.font, visual.style.name)              # set font of node name (Initial Default UNKNOWN (not TNR or Arial))
  setEdgeFontSizeDefault(edge.font.size, visual.style.name)            # set font size of edge (Initial Default 12)
  setNodeFontSizeDefault(node.font.size, visual.style.name)            # set font size of node name (Initial Default 12)
  # shape/style
  setEdgeLineStyleDefault(edge.line.style, visual.style.name)          # "PARALLEL_LINES", "MARQUEE_EQUAL", "DOT", "EQUAL_DASH", "LONG_DASH", "CONTIGUOUS_ARROW", "MARQUEE_DASH", "DASH_DOT", "BACKWARD_SLASH", "FORWARD_SLASH", "VERTICAL_SLASH", "SOLID", "SEPARATE_ARROW", "MARQUEE_DASH_DOT", "ZIGZAG", "SINEWAVE"
  setEdgeSourceArrowShapeDefault(source.arrow, visual.style.name)      # "DELTA", "DIAMOND", "OPEN_CIRCLE", "CIRCLE", "OPEN_HALF_CIRCLE", "CROSS_OPEN_DELTA", "DELTA_SHORT_1", "CROSS_DELTA", "OPEN_DELTA", "OPEN_DIAMOND"
  setEdgeTargetArrowShapeDefault(target.arrow, visual.style.name)      # "DIAMOND_SHORT_1", "DELTA_SHORT_2", "OPEN_SQUARE", "NONE", "SQUARE", "DIAMOND_SHORT_2", "T", "HALF_BOTTOM", "HALF_TOP", "ARROW_SHORT", "HALF_CIRCLE"
  setNodeShapeDefault(node.shape, visual.style.name)                   # "RECTANGLE", "VEE", "OCTAGON", "ELLIPSE", "PARALLELOGRAM", "DIAMOND", "HEXAGON", "ROUND_RECTANGLE", "TRIANGLE"
  # size
  setNodeSizeDefault(node.size, visual.style.name)                     # set size of node; height and width assumed the same
  setEdgeLineWidthDefault(edge.width, visual.style.name)               # set edge line thickness; Initial Default 2
  setNodeBorderWidthDefault(border.width, visual.style.name)           # set border thickness; Initial Default 0
  # opacity
  setEdgeOpacityDefault(edge.opacity, visual.style.name)               # set opacity of edge; 0 - 255 w 0 --> translucent
  setEdgeLabelOpacityDefault(edge.label.opacity, visual.style.name)    # set opacity of edge label; 0 - 255 w 0 --> translucent
  setNodeBorderOpacityDefault(border.opacity, visual.style.name)       # set opacity of border of node; 0 - 255 w 0 --> translucent
  setNodeFillOpacityDefault(node.fill.opacity, visual.style.name)      # set opacity of interior color of node; 0 - 255 w 0 --> translucent
  setNodeLabelOpacityDefault(node.label.opacity, visual.style.name)    # set opacity of name of node; 0 - 255 w 0 --> translucent

  setNodeColorToRatios <- function(plotcol){
    cf <- getTableColumns('node')
    if(!(plotcol %in% getTableColumnNames('node'))){
      print (getTableColumnNames('node'))
      cat("\n","\n","\t", "Which attribute will set node size and color?")
      plotcol <- as.character(readLines(con = stdin(), n = 1))
    }
    limits <- range(cf[, plotcol])
    node.sizes     = c (135, 130, 108, 75, 35, 75, 108, 130, 135)
    #	RATIO is plotted
    #	Blue is negative: Yellow positive, Green in middle
    #
    size.control.points = c (-20, -5, -1.0, 0.0, 1.0, 5.0, 20.0)
    color.control.points = c (-20.0, -10.0, -5.0, -1.0, 0.0, 1.0, 5.0, 10.0, 20.0)
    if(limits[1] < min(size.control.points)) {
      size.control.points = c (limits[1], -15.0, -5.0, 0.0, 5.0, 15.0, 100.0)
      color.control.points = c (limits[1]-1, -10.0, -5.0, -2.25, 0.0, 2.25, 5.0, 10.0, 100.0)
    }
    if(limits[2] > max(size.control.points)) {
      size.control.points = c (limits[1], -15.0, -5.0, 0.0, 5.0, 15.0, limits[2])
      color.control.points = c (limits[1]-1, -10.0, -5.0, -2.25, 0.0, 2.25, 5.0, 10.0, limits[2]+1)
    }
    ratio.colors = c ('#0099FF', '#007FFF','#00BFFF', '#00CCFF', '#00FFFF', '#00EE00', '#FFFF7E', '#FFFF00', '#FFE600', '#FFD700', '#FFCC00')
    setNodeColorMapping (plotcol, color.control.points, ratio.colors, 'c', style.name = visual.style.name)
    setNodeSelectionColorDefault ( "#CC00FF", style.name = visual.style.name)
  }
  setNodeColorToRatios('score')

  setNodeLabelMapping("id", style.name = visual.style.name)

  setNodeBorderColorMapping('node.type', table.column.values = c("kinase", "phosphatase"), colors = c('#FF0000', '#A020F0'), mapping.type = 'd', style.name = visual.style.name)
  setNodeBorderWidthMapping('node.type', table.column.values = c("kinase", "phosphatase"), widths = c(5, 5), mapping.type = 'd', style.name = visual.style.name)

  setVisualStyle(visual.style.name)

}
