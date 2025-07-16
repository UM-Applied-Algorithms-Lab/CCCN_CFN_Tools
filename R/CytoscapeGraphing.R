# NOTE TO SELF: rewrite node color mapping like edge width mapping? to get the ranges based on the data set quartiles
# SECOND NOTE TO SELF: use QUARTILES rather than summary to get better range quartile(x, 0.whatever)

# helper function
NodeAppMap <- function(visual.style.name){                                                         # maps colors based on score
  cf <- getTableColumns('node')                                                                    # gets the table (kind of want to clean this up and just pass in node.table like below)
  limits <- range(cf[, "score"])                                                                   # gets RANGE of scores from table
  node.sizes <- c (135, 130, 108, 75, 35, 75, 108, 130, 135)                                       # set node sizes
  size.control.points = c (-100.0, -15.0, -5.0, 0.0, 5.0, 15.0, 100.0)
  color.control.points = c (-20.0, -10.0, -5.0, -1.0, 0.0, 1.0, 5.0, 10.0, 20.0)                   # Blue is negative, Yellow positive, Green in middle
  if(limits[1] < min(size.control.points)) {
    size.control.points = c (limits[1], -15.0, -5.0, 0.0, 5.0, 15.0, 100.0)
    color.control.points = c (limits[1]-1, -10.0, -5.0, -2.25, 0.0, 2.25, 5.0, 10.0, 100.0)        # if minimum of scores is lower, decrease the lowest (just use -inf?)
  }
  if(limits[2] > max(size.control.points)) {
    size.control.points = c (limits[1], -15.0, -5.0, 0.0, 5.0, 15.0, limits[2])
    color.control.points = c (limits[1]-1, -10.0, -5.0, -2.25, 0.0, 2.25, 5.0, 10.0, limits[2]+1)  # if max of scores is higher, increase the highest (just use inf?)
  }

  node.colors = c('#0099FF', '#007FFF','#00BFFF', '#00CCFF', '#00FFFF', '#00EE00', '#FFFF7E', '#FFFF00', '#FFE600', '#FFD700', '#FFCC00')    # colors line up to the ranges

  setNodeColorMapping(plotcol, color.control.points, node.colors, 'c', style.name = visual.style.name)                                       # set node color mapping
  lockNodeDimensions('TRUE')                                                                                                                 # width and height locked together
  setNodeSizeMapping (names(cf[plotcol]), size.control.points, node.sizes, 'c')                                                              # set node size mapping
  setNodeSelectionColorDefault ( "#CC00FF", style.name = visual.style.name)                                                                  # set selection color
}





# helper function
 EdgeAppMap <- function(edge.table, visual.style.name){

  edgevalues <- getTableColumns('edge',c('shared name', 'shared interaction', 'weight', 'interaction')) # get edge values
  edgevalues['weight'] <- abs(as.numeric(edgevalues['weight'][[1]]))                                    # make weights pos
  edgevalues['weight'] <- lapply(edgevalues['weight'], function(x) log2(x * 10) + 2)                    # make scale logarithmic

  colnames(edgevalues) <- c('shared name', 'shared interaction', 'width', 'interaction')                # rename it width bc we're doing passthrough mapping

  setEdgeSelectionColorDefault ("#FF69B4")                                                              # selection edge color

  # list of all edge types (you only have SOME in data set, Maddie)
  edgeTypes <- c("PHOSPHORYLATION", "pp", "controls-phosphorylation-of", "controls-expression-of", "controls-transport-of", "controls-state-change-of", "ACETYLATION", "Physical interactions", "BioPlex", "in-complex-with", 'experiments', 'experiments_transferred', 'database', 'database_transfered', "Pathway", "Predicted", "Genetic interactions", "correlation", "negative correlation", "positive correlation", 'combined_score', "merged" , "intersect", "peptide", 'homology', "Shared protein domains")

  # list of colors corresponding to edge types
  edgecolors <- c("#FF0000", "#FF0000", "#FF0000", "#FF00FF", "#EE82EE", "#A020F0", "#FF7F00", "#00FF00", "#00EE00", "#00CD00", "#76EEC6", "#76EEC6", "#00FFFF", "#00FFFF", "#00E5EE", "#00EEEE", "#20B2AA", "#FFD700", "#0000FF", "#FFFF00", "#708090", "#2F4F4F", "#BEBEBE", "#000000", "#FFA500", "#EE9A00")

  # list of arrow options corresponding to edge types
  myarrows <- c ('Arrow', 'Arrow', 'Arrow', 'Arrow', 'Arrow', 'Arrow', "Arrow", 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None')

  # Fix phosphorylation edge color to red with merged edges  < -- not sure if these comments are still relevant
  # work around for misaligned mapping                       < -- ignoring for now Ig

  # OH if it's in shared interaction, make the interaction column that
  # shouldn't this already be so?

  if (length(edgevalues[grep("pp", edgevalues$'shared interaction'), 'interaction']) > 0) {                       # make interaction of pp shared interaction "PHOSPHORYLATION"
    edgevalues[grep("pp", edgevalues$'shared interaction'), 'interaction'] <- "PHOSPHORYLATION"}

  if (length(edgevalues[grep("PHOSPHORYLATION", edgevalues$'shared interaction'), 'interaction']) > 0) {          # make interaction of PHOSPHORYLATION shared interaction "PHOSPHORYLATION"
    edgevalues[grep("PHOSPHORYLATION", edgevalues$'shared interaction'), 'interaction'] <- "PHOSPHORYLATION"}

  if (length(edgevalues[grep("phosphorylation", edgevalues$'shared interaction'), 'interaction']) > 0) {          # make interaction of phosphorylation shared interaction "PHOSPHORYLATION"
    edgevalues[grep("phosphorylation", edgevalues$'shared interaction'), 'interaction'] <- "PHOSPHORYLATION"}

  if (length(edgevalues[grep("ACETYLATION", edgevalues$'shared interaction'), 'interaction']) > 0) {              # make interaction of ACETYLATION shared interaction "ACETYLATION"
    edgevalues[grep("ACETYLATION", edgevalues$'shared interaction'), 'interaction'] <- "ACETYLATION"}

  loadTableData(edgevalues, table = 'edge', table.key.column = 'SUID')                                            # load the table

  matchArrowColorToEdge('TRUE')                                                                                   # match arrow color to edge color (love RCy3 for this)

  setEdgeLineWidthMapping('width', mapping.type = 'p', style.name = visual.style.name)                            # map the edge width
  setEdgeTargetArrowMapping('interaction', edgeTypes, myarrows, default.shape='None')                             # map the target arrow
  setEdgeColorMapping('interaction', edgeTypes, edgecolors, 'd', default.color="#FFFFFF")                         # map the edge color

}





# helper function
NodeBorderMapping <- function(visual.style.name){

  # list of molecule classes
  molclasses <- c("unknown", "receptor tyrosine kinase",  "SH2 protein", "SH2-SH3 protein", "SH3 protein", "tyrosine kinase",  "SRC-family kinase", "kinase", "phosphatase", "transcription factor", "RNA binding protein")
  # list of node shapes
  nodeshapes <- c("ELLIPSE","ROUND_RECTANGLE", "VEE", "VEE", "TRIANGLE", "HEXAGON", "DIAMOND", "OCTAGON", "OCTAGON", "PARALLELOGRAM", "RECTANGLE")

  # map molecule classes to node shapes
  setNodeShapeMapping ("node.type", molclasses, nodeshapes, default.shape="ELLIPSE", style.name = visual.style.name)
  # map node.type with this list of molecule types to different widths
  setNodeBorderWidthMapping("node.type", c("acetyltransferase", "methyltransferase", "membrane protein", "receptor tyrosine kinase", "G protein-coupled receptor", "SRC-family kinase", "tyrosine kinase", "kinase", "phosphatase"),
                            widths=c(12,12,8,16,16,12,12,12,14), mapping.type = 'd', default.width = 4, style.name = visual.style.name)

  # map the different molecule types to these different node border colors
  setNodeBorderColorMapping(
    table.column = "node.type",
    table.column.values = c("deacetylase", "acetyltransferase", "demethylase", "methyltransferase", "membrane protein", "kinase", "tyrosine kinase", "SRC-family kinase", "phosphatase", "tyrosine phosphatase", "G protein-coupled receptor", "receptor tirosine kinase"),
    colors = c("#FF8C00", "#FF8C00", "#005CE6", "#005CE6", "#6600CC", "#EE0000", "#EE0000", "#EE0000", "#FFEC8B", "#FFEC8B", "#BF3EFF", "#BF3EFF"),
    mapping.type = 'd',
    default.color = '#abb2b9',
    style.name = visual.style.name
  )
}





# helper function
SetStandards <- function(visual.style.name, background.color, edge.label.color, edge.line.color, node.border.color, node.label.color,
                         default.font, node.font.size, edge.font.size, edge.line.style, source.arrow, target.arrow, node.size, edge.width, border.width,
                         edge.opacity, edge.label.opacity, border.opacity, node.label.opacity, node.fill.opacity
                         ){

  setEdgeSourceArrowColorDefault(edge.line.color, visual.style.name)   # set color of arrow from source
  setEdgeTargetArrowColorDefault(edge.line.color, visual.style.name)   # set color of arrow to target
  setNodeLabelPositionDefault("C", "C", "c", 0, 0, visual.style.name)  # What part of the node label is aligned "C", "NW", "N", "NE", "E", "SE", "S", "SW", "W"
                                                                       # to what part of the node graphic       "C", "NW", "N", "NE", "E", "SE", "S", "SW", "W"
                                                                       # "l", "r", "c"
                                                                       # amount offset in the x direction
                                                                       # amount offset in the y direction

  # CUSTOMIZED: user inputs can change the following:

  # colors
  setBackgroundColorDefault(background.color, visual.style.name)       # set color of background
  setEdgeLabelColorDefault(edge.label.color, visual.style.name)        # set color of edge label
  setEdgeColorDefault(edge.line.color, visual.style.name)              # set color of edge
  setNodeLabelColorDefault(node.label.color, visual.style.name)        # set color of node name
  # fonts
  setEdgeFontFaceDefault(default.font, visual.style.name)              # set font of edge
  setNodeFontFaceDefault(default.font, visual.style.name)              # set font of node name (Initial Default UNKNOWN (not TNR or Arial))
  setEdgeFontSizeDefault(edge.font.size, visual.style.name)            # set font size of edge (Initial Default 12)
  setNodeFontSizeDefault(node.font.size, visual.style.name)            # set font size of node name (Initial Default 12)
  # shape/style
  setEdgeLineStyleDefault(edge.line.style, visual.style.name)          # "PARALLEL_LINES", "MARQUEE_EQUAL", "DOT", "EQUAL_DASH", "LONG_DASH", "CONTIGUOUS_ARROW", "MARQUEE_DASH", "DASH_DOT", "BACKWARD_SLASH", "FORWARD_SLASH", "VERTICAL_SLASH", "SOLID", "SEPARATE_ARROW", "MARQUEE_DASH_DOT", "ZIGZAG", "SINEWAVE"
  setEdgeSourceArrowShapeDefault(source.arrow, visual.style.name)      # "DELTA", "DIAMOND", "OPEN_CIRCLE", "CIRCLE", "OPEN_HALF_CIRCLE", "CROSS_OPEN_DELTA", "DELTA_SHORT_1", "CROSS_DELTA", "OPEN_DELTA", "OPEN_DIAMOND"
  setEdgeTargetArrowShapeDefault(target.arrow, visual.style.name)      # "DIAMOND_SHORT_1", "DELTA_SHORT_2", "OPEN_SQUARE", "NONE", "SQUARE", "DIAMOND_SHORT_2", "T", "HALF_BOTTOM", "HALF_TOP", "ARROW_SHORT", "HALF_CIRCLE"
  # size
  setNodeSizeDefault(node.size, visual.style.name)                     # set size of node; height and width assumed the same
  # opacity
  setEdgeOpacityDefault(edge.opacity, visual.style.name)               # set opacity of edge; 0 - 255 w 0 --> translucent
  setEdgeLabelOpacityDefault(edge.label.opacity, visual.style.name)    # set opacity of edge label; 0 - 255 w 0 --> translucent
  setNodeBorderOpacityDefault(border.opacity, visual.style.name)       # set opacity of border of node; 0 - 255 w 0 --> translucent
  setNodeFillOpacityDefault(node.fill.opacity, visual.style.name)      # set opacity of interior color of node; 0 - 255 w 0 --> translucent
  setNodeLabelOpacityDefault(node.label.opacity, visual.style.name)    # set opacity of name of node; 0 - 255 w 0 --> translucent

}





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
#'
#' @param default.font Font style of edge and node names; defaults to "Times New Roman"
#' @param node.font.size Font size of the node name; defaults to 12
#' @param edge.font.size Font size of the edge name; defaults to 8
#'
#' @param edge.line.style Type of edge style; defaults to "SOLID"; options include: "PARALLEL_LINES", "MARQUEE_EQUAL", "DOT", "EQUAL_DASH", "LONG_DASH", "CONTIGUOUS_ARROW", "MARQUEE_DASH", "DASH_DOT", "BACKWARD_SLASH", "FORWARD_SLASH", "VERTICAL_SLASH", "SOLID", "SEPARATE_ARROW", "MARQUEE_DASH_DOT", "ZIGZAG", "SINEWAVE"
#' @param source.arrow Type of arrow coming from the source gene; defaults to "NONE"; options include: "DELTA", "DIAMOND", "OPEN_CIRCLE", "CIRCLE", "OPEN_HALF_CIRCLE", "CROSS_OPEN_DELTA", "DELTA_SHORT_1", "CROSS_DELTA", "OPEN_DELTA", "OPEN_DIAMOND", "DIAMOND_SHORT_1", "DELTA_SHORT_2", "OPEN_SQUARE", "NONE", "SQUARE", "DIAMOND_SHORT_2", "T", "HALF_BOTTOM", "HALF_TOP", "ARROW_SHORT", "HALF_CIRCLE"
#' @param target.arrow Type of arrow going to the target gene; defaults to "NONE"; options include: "DELTA", "DIAMOND", "OPEN_CIRCLE", "CIRCLE", "OPEN_HALF_CIRCLE", "CROSS_OPEN_DELTA", "DELTA_SHORT_1", "CROSS_DELTA", "OPEN_DELTA", "OPEN_DIAMOND", "DIAMOND_SHORT_1", "DELTA_SHORT_2", "OPEN_SQUARE", "NONE", "SQUARE", "DIAMOND_SHORT_2", "T", "HALF_BOTTOM", "HALF_TOP", "ARROW_SHORT", "HALF_CIRCLE"
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
                     background.color = '#faf1dd', edge.label.color = '#17202a', edge.line.color = '#abb2b9', node.border.color = '#145a32', node.label.color = '#000000',
                     default.font = "Times New Roman", node.font.size = 12, edge.font.size = 8,
                     edge.line.style = 'SOLID', source.arrow = 'NONE', target.arrow = 'NONE',
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

  genes <- unique(c(cfn$Gene.1, cfn$Gene.2))                                           # get unique names of genes

  ptmnew <- ptmtable[, c("PTM", "HCC4006_Erlotinib")]                                  # just take ptms and another col doesn't matter bc we replace
  ptmnew$PTM <- sapply(ptmtable$PTM, function(x) strsplit(x, split = " ")[[1]][1])     # get just the gene names in the PTM column
  colnames(ptmnew) <- c("Gene", "score")                                               # rename cols

  ptmnew$score <- sapply(1:length(rownames(ptmtable)), function(i) min(as.numeric(ptmtable[-1][i, ]), na.rm = TRUE))    # TAKE MINIMUM SCORE ACROSS PTM COLS

  cfn.edges <- data.frame(matrix(data = 0, nrow = length(rownames(cfn)), ncol = 4), stringsAsFactors = FALSE)           # initialize empty edge df
  cfn.nodes <- data.frame(matrix(data = 0, nrow = length(genes), ncol = 3), stringsAsFactors = FALSE)                   # initialize empty node df

  colnames(cfn.edges) <- c("source", "target", "interaction", "weight")                                                 # name cols of edge and node tables
  colnames(cfn.nodes) <- c("id", "node.type", "score")

  cfn.edges$source <- cfn$Gene.1                                                            # enter vals for edge table
  cfn.edges$target <- cfn$Gene.2
  cfn.edges$interaction <- cfn$Interaction
  cfn.edges$weight <- format(as.numeric(cfn$PPI.weight), scientific = FALSE, trim = TRUE)   # NOTE: formatting is very important for the edge line width

  # enter vals for node table
  cfn.nodes$id <- genes
  cfn.nodes$node.type <- as.character(sapply(cfn.nodes$id, function(x) funckey$nodeType[which(funckey$Gene.Name == x)]))   # steal node type from funckey
  cfn.nodes$score <- as.numeric(sapply(cfn.nodes$id, function(x) sum(ptmnew$score[which(ptmnew$Gene == x)])))              # steal score from ptmnew

  cyscape <- createNetworkFromDataFrames(cfn.nodes, cfn.edges, title = Network.title, collection = Network.collection)     # create network (not sure if storing it does anything?)

  copyVisualStyle("default", visual.style.name)   # create visual style



  # CUSTOMIZATION FROM HERE ON OUT

  setNodeLabelMapping("id", style.name = visual.style.name)   # make the label names appear

  NodeAppMap(visual.style.name)                         # map the node colors CURRENTLY NOT SIZES

  EdgeAppMap(cfn.edges, visual.style.name)              # map edge width CURRENTLY NOT COLORS

  NodeBorderMapping(visual.style.name)                        # map node border (color, thickness)

  # we got standards in this joint
  SetStandards(visual.style.name, background.color, edge.label.color, edge.line.color, node.border.color, node.label.color,
               default.font, node.font.size, edge.font.size, edge.line.style, source.arrow, target.arrow, node.size, edge.width, border.width,
               edge.opacity, edge.label.opacity, border.opacity, node.label.opacity, node.fill.opacity)

  setVisualStyle(visual.style.name)                           # set vis style

}



# Graph CoCluster Correlation Network
#
# Creates a cytoscape graph of the cocluster correlation network. Ensure that you have the cytoscape app open and the RCy3 package downloaded and libraried.
#
# If error occurs, run cytoscapePing() to ensure connection to the Cytoscape Interface and try again.
#
# @param ptm.cccn A matrix showing strength of relationships between post-translational modifications using common clusters between the three distance metrics
# @param ptmtable A dataset for post-translational modifications. Formatted with numbered rows, and the first column containing PTM names. The rest of the column names should be drugs. Values are numeric values that represent how much the PTM has reacted to the drug.
# @param funckey A table graphing gene names to type of protein; defaults to internal database at cccn.cfn.tools::ex.funckey
# @param Network.title Desired title for the created Cytoscape Network; defaults to "cfn"
# @param Network.collection Desired name for the collection created on Cytoscape in which the network will reside; defaults to "cccn.cfn.tools"
# @param visual.style.name Desired name for the visual style created on Cytoscape; defaults to "cccn.cfn.tools.style"
#
# @param background.color Hex code of background color of graph; defaults to '#fcf3cf'
# @param edge.label.color Hex code of edge label color of graph; defaults to '#17202a'
# @param edge.line.color Hex code of edge line color of graph; defaults to '#abb2b9'
# @param node.border.color Hex code of node border color of graph; defaults to '#145a32'
# @param node.label.color Hex code of node label color of graph; defaults to '#145a32'
#
# @param default.font Font style of edge and node names; defaults to "Times New Roman"
# @param node.font.size Font size of the node name; defaults to 12
# @param edge.font.size Font size of the edge name; defaults to 8
#
# @param edge.line.style Type of edge style; defaults to "SOLID"; options include: "PARALLEL_LINES", "MARQUEE_EQUAL", "DOT", "EQUAL_DASH", "LONG_DASH", "CONTIGUOUS_ARROW", "MARQUEE_DASH", "DASH_DOT", "BACKWARD_SLASH", "FORWARD_SLASH", "VERTICAL_SLASH", "SOLID", "SEPARATE_ARROW", "MARQUEE_DASH_DOT", "ZIGZAG", "SINEWAVE"
# @param source.arrow Type of arrow coming from the source gene; defaults to "NONE"; options include: "DELTA", "DIAMOND", "OPEN_CIRCLE", "CIRCLE", "OPEN_HALF_CIRCLE", "CROSS_OPEN_DELTA", "DELTA_SHORT_1", "CROSS_DELTA", "OPEN_DELTA", "OPEN_DIAMOND", "DIAMOND_SHORT_1", "DELTA_SHORT_2", "OPEN_SQUARE", "NONE", "SQUARE", "DIAMOND_SHORT_2", "T", "HALF_BOTTOM", "HALF_TOP", "ARROW_SHORT", "HALF_CIRCLE"
# @param target.arrow Type of arrow going to the target gene; defaults to "NONE"; options include: "DELTA", "DIAMOND", "OPEN_CIRCLE", "CIRCLE", "OPEN_HALF_CIRCLE", "CROSS_OPEN_DELTA", "DELTA_SHORT_1", "CROSS_DELTA", "OPEN_DELTA", "OPEN_DIAMOND", "DIAMOND_SHORT_1", "DELTA_SHORT_2", "OPEN_SQUARE", "NONE", "SQUARE", "DIAMOND_SHORT_2", "T", "HALF_BOTTOM", "HALF_TOP", "ARROW_SHORT", "HALF_CIRCLE"
#
# @param node.size Size of the node; defaults to 50. PLEASE NOTE: width and height can be changed independently using RCy3 directly (first run lockNodeDimensions(FALSE) and then setNodeWidthDefault() and setNodeHeightDefault())
# @param edge.width Width of the edge line; defaults to 2
# @param border.width Width of the border of the node; defaults to 1
#
# @param edge.opacity Opacity of the edge line on a scale of 0 - 255 with 0 being transparent; defaults to 175
# @param edge.label.opacity Opacity of the edge label on a scale of 0 - 255 with 0 being transparent; defaults to 255
# @param border.opacity Opacity of the node border on a scale of 0 - 255 with 0 being transparent; defaults to 255
# @param node.label.opacity Opacity of the node label on a scale of 0 - 255 with 0 being transparent; defaults to 255
# @param node.fill.opacity Opacity of the node fill on a scale of 0 - 255 with 0 being transparent; defaults to 255
#
# @return A cytoscape graph of the cluster filtered network
# @export
#
# @examples
# # GraphCFN(ex.cfn)
# # See vignette for default graph
#GraphPTMCccn <- function(ptm.cccn, ptmtable, funckey = cccn.cfn.tools::ex.funckey, Network.title = "PTM.cccn", Network.collection = "cccn.cfn.tools", visual.style.name = "cccn.cfn.tools.style",
#                     background.color = '#faf1dd', edge.label.color = '#17202a', edge.line.color = '#abb2b9', node.border.color = '#145a32', node.label.color = '#000000',
#                     default.font = "Times New Roman", node.font.size = 12, edge.font.size = 8,
#                     edge.line.style = 'SOLID', source.arrow = 'NONE', target.arrow = 'NONE',
#                     node.size = 50, edge.width = 2, border.width = 1,
#                     edge.opacity = 175, edge.label.opacity = 255, border.opacity = 255, node.label.opacity = 255, node.fill.opacity = 255
#                     ){
#
#  if(!(exists("RCy3"))){                      # check if RCy3 is libraried
#    if(system.file(package="RCy3") == ""){    # check if RCy3 is installed at all
#      BiocManager::install(version='devel')   # The following initializes usage of Bioc devel
#      BiocManager::install("RCy3")            # install!
#    }
#    library("RCy3")                           # library
#    cytoscapePing()                           # ensure connection?
#  }
#
#  # ACTUAL CODE AND DATA PROCESSING
#
#  ptms.list <- unique(c(rownames(ptm.cccn), colnames(ptm.cccn)))                       # get unique names of ptms
#  ptm.cccn[which(is.na(ptm.cccn))] <- 0                                                # remove NAs
#
#  ptmnew <- ptmtable[, c("PTM", "HCC4006_Erlotinib")]                                  # just take ptms and another col doesn't matter bc we replace
#  colnames(ptmnew) <- c("PTM", "score")                                                # rename cols
#
#  ptmnew$score <- sapply(1:length(rownames(ptmtable)), function(i) min(as.numeric(ptmtable[-1][i, ]), na.rm = TRUE))         # TAKE MINIMUM SCORE ACROSS PTM COLS
#
#  ptm.cccn.edges <- data.frame(matrix(data = 0, nrow = length(rownames(ptm.cccn))**2, ncol = 4), stringsAsFactors = FALSE)   # initialize empty edge df
#  ptm.cccn.nodes <- data.frame(matrix(data = 0, nrow = length(ptms.list), ncol = 3), stringsAsFactors = FALSE)                # initialize empty node df
#
#  colnames(ptm.cccn.edges) <- c("source", "target", "interaction", "weight")                                                 # name cols of edge and node tables
#  colnames(ptm.cccn.nodes) <- c("id", "node.type", "score")
#
#  edge.grid <- expand.grid(
#    row = rownames(ptm.cccn),
#    col = colnames(ptm.cccn),
#    stringsAsFactors = FALSE
#  )
#  values <- mapply(function(r, c) ptm.cccn[r, c], edge.grid$row, edge.grid$col)
#
#  ptm.cccn.edges$source <- edge.grid[, 1]
#  ptm.cccn.edges$target <- edge.grid[, 2]
#  ptm.cccn.edges$interaction <- sapply(values, function(x) if(x > 0){return("positive")} else {return("negative")})
#  ptm.cccn.edges$weight <- sapply(values, function(x)abs(x))
#
#  ptm.cccn.edges <- ptm.cccn.edges[which(!(ptm.cccn.edges$weight == 0)), ]
#  ptm.cccn.edges$weight <- format(as.numeric(ptm.cccn.edges$weight), scientific = FALSE, trim = TRUE)   # NOTE: formatting is very important for the edge line width
#
#  # enter vals for node table
#  ptm.cccn.nodes$id <- ptms.list
#  ptm.cccn.nodes$node.type <- as.character(sapply(ptm.cccn.nodes$id, function(x) funckey$nodeType[which(funckey$Gene.Name == strsplit(x, split = ' ')[[1]][1])]))   # steal node type from funckey
#  ptm.cccn.nodes$score <- as.numeric(sapply(ptm.cccn.nodes$id, function(x) sum(ptmnew$score[which(ptmnew$Gene == strsplit(x, split = ' ')[[1]][1])])))              # steal score from ptmnew
#
#  cyscape <- createNetworkFromDataFrames(ptm.cccn.nodes, ptm.cccn.edges, title = Network.title, collection = Network.collection)     # create network (not sure if storing it does anything?)
#
#  copyVisualStyle("default", visual.style.name)   # create visual style
#
#
#
#  # CUSTOMIZATION FROM HERE ON OUT
#
#  setNodeLabelMapping("id", style.name = visual.style.name)   # make the label names appear
#
#  NodeAppMap('score', visual.style.name)                # map the node colors CURRENTLY NOT SIZES
#
#   EdgeAppMap(ptm.cccn.edges, visual.style.name)              # map edge width CURRENTLY NOT COLORS
#
#  NodeBorderMapping(visual.style.name)                        # map node border (color, thickness)

#  # we got standards in this joint
#  SetStandards(visual.style.name, background.color, edge.label.color, edge.line.color, node.border.color, node.label.color,
#               default.font, node.font.size, edge.font.size, edge.line.style, source.arrow, target.arrow, node.size, edge.width, border.width,
#               edge.opacity, edge.label.opacity, border.opacity, node.label.opacity, node.fill.opacity)
#
#  setVisualStyle(visual.style.name)                           # set vis style
#
#}
