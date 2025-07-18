# helper function
NodeAppMap <- function(visual.style.name){                                                         # maps colors based on score
  cf <- RCy3::getTableColumns('node')                                                              # gets the table
  limits <- range(cf[, "score"])                                                                   # gets RANGE of scores from table
  node.sizes <- c(135, 130, 108, 75, 35, 75, 108, 130, 135)                                        # set node sizes
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

  node.colors = c('#0099FF', '#007FFF','#00BFFF', '#00CCFF', '#00FFFF', '#00EE00', '#FFFF7E', '#FFFF00', '#FFE600', '#FFD700', '#FFCC00')  # colors line up to the ranges

  RCy3::setNodeColorMapping("score", color.control.points, node.colors, 'c', style.name = visual.style.name)                               # set node color mapping
  RCy3::setNodeSelectionColorDefault ( "#CC00FF", style.name = visual.style.name)                                                          # set selection color
  RCy3::lockNodeDimensions('TRUE', style = visual.style.name)                                                                              # width and height locked together

  RCy3::setNodeSizeMapping(table.column = "score", table.column.values = size.control.points, sizes = node.sizes, mapping.type = 'c', default.size = 103, style.name = visual.style.name) # set node size mapping
}





# helper function
 EdgeAppMap <- function(visual.style.name){

  edgevalues <- RCy3::getTableColumns('edge',c('shared name', 'shared interaction', 'weight', 'interaction')) # get edge values
  edgevalues['weight'] <- abs(as.numeric(edgevalues['weight'][[1]]))                                          # make weights pos
  edgevalues['weight'] <- lapply(edgevalues['weight'], function(x) log2(x * 10) + 2)                          # make scale logarithmic

  colnames(edgevalues) <- c('shared name', 'shared interaction', 'width', 'interaction')                      # rename it width bc we're doing passthrough mapping

  RCy3::setEdgeSelectionColorDefault ("#FF69B4", style.name = visual.style.name)                              # selection edge color

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

  RCy3::loadTableData(edgevalues, table = 'edge', table.key.column = 'SUID')                                      # load the table

  RCy3::matchArrowColorToEdge('TRUE', style.name = visual.style.name)                                             # match arrow color to edge color (love RCy3 for this)

  RCy3::setEdgeLineWidthMapping('width', mapping.type = 'p', style.name = visual.style.name)                                     # map the edge width
  RCy3::setEdgeTargetArrowMapping('interaction', edgeTypes, myarrows, default.shape='None', style.name = visual.style.name)      # map the target arrow
  RCy3::setEdgeColorMapping('interaction', edgeTypes, edgecolors, 'd', default.color="#FFFFFF", style.name = visual.style.name)  # map the edge color

}





# helper function
NodeBorderMapping <- function(visual.style.name){

  # list of molecule classes
  molclasses <- c("unknown", "receptor tyrosine kinase",  "SH2 protein", "SH2-SH3 protein", "SH3 protein", "tyrosine kinase",  "SRC-family kinase", "kinase", "phosphatase", "transcription factor", "RNA binding protein")
  # list of node shapes
  nodeshapes <- c("ELLIPSE","ROUND_RECTANGLE", "VEE", "VEE", "TRIANGLE", "HEXAGON", "DIAMOND", "OCTAGON", "OCTAGON", "PARALLELOGRAM", "RECTANGLE")

  # map molecule classes to node shapes
  RCy3::setNodeShapeMapping ("node.type", molclasses, nodeshapes, default.shape="ELLIPSE", style.name = visual.style.name)
  # map node.type with this list of molecule types to different widths
  RCy3::setNodeBorderWidthMapping("node.type", c("acetyltransferase", "methyltransferase", "membrane protein", "receptor tyrosine kinase", "G protein-coupled receptor", "SRC-family kinase", "tyrosine kinase", "kinase", "phosphatase"),
                            widths=c(12,12,8,16,16,12,12,12,14), mapping.type = 'd', default.width = 4, style.name = visual.style.name)

  # map the different molecule types to these different node border colors
  RCy3::setNodeBorderColorMapping(
    table.column = "node.type",
    table.column.values = c("deacetylase", "acetyltransferase", "demethylase", "methyltransferase", "membrane protein", "kinase", "tyrosine kinase", "SRC-family kinase", "phosphatase", "tyrosine phosphatase", "G protein-coupled receptor", "receptor tyrosine kinase"),
    colors = c("#FF8C00", "#FF8C00", "#005CE6", "#005CE6", "#6600CC", "#EE0000", "#EE0000", "#EE0000", "#FFEC8B", "#FFEC8B", "#BF3EFF", "#BF3EFF"),
    mapping.type = 'd',
    default.color = '#abb2b9',
    style.name = visual.style.name
  )
}





# helper function
SetStandards <- function(visual.style.name,
                         background.color, edge.label.color, node.label.color, default.font, node.font.size, edge.font.size,
                         edge.line.style,edge.opacity, edge.label.opacity, border.opacity, node.label.opacity, node.fill.opacity
                         ){

  RCy3::setNodeLabelPositionDefault("C", "C", "c", 0, 0, visual.style.name)  # What part of the node label is aligned "C", "NW", "N", "NE", "E", "SE", "S", "SW", "W"
                                                                             # to what part of the node graphic       "C", "NW", "N", "NE", "E", "SE", "S", "SW", "W"
                                                                             # "l", "r", "c"
                                                                             # amount offset in the x direction
                                                                             # amount offset in the y direction

  # CUSTOMIZED: user inputs can change the following:

  # colors
  RCy3::setBackgroundColorDefault(background.color, visual.style.name)       # set color of background
  RCy3::setEdgeLabelColorDefault(edge.label.color, visual.style.name)        # set color of edge label
  RCy3::setNodeLabelColorDefault(node.label.color, visual.style.name)        # set color of node name
  # fonts
  RCy3::setEdgeFontFaceDefault(default.font, visual.style.name)              # set font of edge
  RCy3::setNodeFontFaceDefault(default.font, visual.style.name)              # set font of node name (Initial Default UNKNOWN (not TNR or Arial))
  RCy3::setEdgeFontSizeDefault(edge.font.size, visual.style.name)            # set font size of edge (Initial Default 12)
  RCy3::setNodeFontSizeDefault(node.font.size, visual.style.name)            # set font size of node name (Initial Default 12)
  # shape/style
  RCy3::setEdgeLineStyleDefault(edge.line.style, visual.style.name)          # "PARALLEL_LINES", "MARQUEE_EQUAL", "DOT", "EQUAL_DASH", "LONG_DASH", "CONTIGUOUS_ARROW", "MARQUEE_DASH", "DASH_DOT", "BACKWARD_SLASH", "FORWARD_SLASH", "VERTICAL_SLASH", "SOLID", "SEPARATE_ARROW", "MARQUEE_DASH_DOT", "ZIGZAG", "SINEWAVE"
  # opacity
  RCy3::setEdgeOpacityDefault(edge.opacity, visual.style.name)               # set opacity of edge; 0 - 255 w 0 --> translucent
  RCy3::setEdgeLabelOpacityDefault(edge.label.opacity, visual.style.name)    # set opacity of edge label; 0 - 255 w 0 --> translucent
  RCy3::setNodeBorderOpacityDefault(border.opacity, visual.style.name)       # set opacity of border of node; 0 - 255 w 0 --> translucent
  RCy3::setNodeFillOpacityDefault(node.fill.opacity, visual.style.name)      # set opacity of interior color of node; 0 - 255 w 0 --> translucent
  RCy3::setNodeLabelOpacityDefault(node.label.opacity, visual.style.name)    # set opacity of name of node; 0 - 255 w 0 --> translucent

}





#' @title Graph Cluster Filtered Network
#'
#' @description Creates a cytoscape graph of the cluster filtered network. Ensure that you have the Cytoscape app open and the RCy3 package downloaded.
#' If RCy3 is installed and the Cytoscape App is open yet an error still occurs, run cytoscapePing() to ensure connection to the Cytoscape Interface and try again.
#' The package RCy3 is required for this function. To download, run:
#' - if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#' - BiocManager::install("RCy3")
#'
#' @details Graph Key
#' - Node Size
#'   - Greater the node size, larger the absolute value of the score
#' - Blue Node
#'   - Negative score
#' - Yellow Node
#'   - Positive score
#' - Green Node
#'   - Approximately zero score
#' Node Shapes:
#' - "ELLIPSE"
#'   - unknown
#' - "ROUND_RECTANGLE"
#'   - receptor tyrosine kinase
#' - "VEE"
#'   - SH2 protein
#'   or
#'   - SH2-SH3 protein
#' - "TRIANGLE"
#'   - SH3 protein
#' - "HEXAGON"
#'   - tyrosine kinase
#' - "DIAMOND"
#'   - SRC-family kinase
#' - "OCTAGON"
#'   - kinase
#'   or
#'   - phosphatase
#' - "PARALLELOGRAM"
#'   - transcription factor
#' - "RECTANGLE"
#'   - RNA binding protein
#' Node Border Colors:
#' - Orange
#'   - deacetylase
#'   or
#'   - acetyltransferase
#' - Blue
#'   - demethylase
#'   or
#'   - methyltransferase
#' - Royal Purple
#'   - membrane protein
#' - Red
#'   - kinase
#'   or
#'   - tyrosine kinase
#'   or
#'   - SRC-family kinase
#' - Yellow
#'   - phosphatase
#'   or
#'   - tyrosine phosphatase
#' - Lilac
#'   - G protein-coupled receptor
#'   or
#'   - receptor tyrosine kinase
#' - Grey
#'   - default
#' - Edge Thickness
#'   - Stronger correlation
#'Edge Colors:
#' - Red
#'   - Phosphorylation
#'   or
#'   - pp
#'   or
#'   - controls-phosphorylation-of
#' - Bright Magenta
#'   - controls-expression-of
#' - Dull Magenta
#'   - controls-transport-of
#' - Purple
#'   - controls-state-change-of
#' - Blood Orange
#'   - Acetylation
#' - Lime Green
#'   - Phystical interactions
#' - Green
#'   - BioPlex
#' - Dull Green
#'   - in-complex-with
#' - Seafoam Green
#'   - experiments
#'   or
#'   - experiments_transferred
#' - Cyan
#'   - database
#'   or
#'   - database_transferred
#' - Teal
#'   - Pathway
#'   or
#'   - Predicted
#' - Dark Turquoise
#'   - Genetic interactions
#' - Yellow-Orange
#'   - correlation
#' - Royal Blue
#'   - negative correlation
#' - Bright Yellow
#'   - positive correlation
#' - Grey
#'   - combined_score
#' - Dark Grey
#'   - merged
#' - Light Grey
#'   - intersect
#' - Black
#'   - peptide
#' - Orange
#'   - homology
#' - Dull Orange
#'   - Shared protein domains
#' - White
#'   - Default
#' Arrow Types:
#' - Arrow
#'   - Phosphorylation
#'   or
#'   - pp
#'   or
#'   - controls-phosphorylation-of
#'   or
#'   - controls-expression-of
#'   or
#'   - controls-transport-of
#'   or
#'   - controls-state-change-of
#'   or
#'   - Acetylation
#' - No Arrow
#'   - Default
#'
#' @param cfn A version of ppi.network with only the edges that exist in cccn.matrix and have non-zero weights
#' @param ptmtable A dataset for post-translational modifications. Formatted with numbered rows, and the first column containing PTM names. The rest of the column names should be drugs. Values are numeric values that represent how much the PTM has reacted to the drug.
#' @param funckey A table graphing gene names to type of protein; defaults to internal database at PTMsToPathways::ex.funckey
#' @param Network.title Desired title for the created Cytoscape Network; defaults to "cfn"
#' @param Network.collection Desired name for the collection created on Cytoscape in which the network will reside; defaults to "PTMsToPathways"
#' @param visual.style.name Desired name for the visual style created on Cytoscape; defaults to "PTMsToPathways.style"
#'
#' @param background.color Hex code of background color of graph; defaults to '#fcf3cf'
#' @param edge.label.color Hex code of edge label color of graph; defaults to '#17202a'
#' @param node.label.color Hex code of node label color of graph; defaults to '#145a32'
#'
#' @param default.font Font style of edge and node names; defaults to "Times New Roman"
#' @param node.font.size Font size of the node name; defaults to 12
#' @param edge.font.size Font size of the edge name; defaults to 8
#'
#' @param edge.line.style Type of edge style; defaults to "SOLID"; options include: "PARALLEL_LINES", "MARQUEE_EQUAL", "DOT", "EQUAL_DASH", "LONG_DASH", "CONTIGUOUS_ARROW", "MARQUEE_DASH", "DASH_DOT", "BACKWARD_SLASH", "FORWARD_SLASH", "VERTICAL_SLASH", "SOLID", "SEPARATE_ARROW", "MARQUEE_DASH_DOT", "ZIGZAG", "SINEWAVE"
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
GraphCfn <- function(cfn, ptmtable, funckey = PTMsToPathways::ex.funckey, Network.title = "cfn", Network.collection = "PTMsToPathways", visual.style.name = "PTMsToPathways.style",
                     background.color = "#949494", edge.label.color = '#17202a', node.label.color = '#000000',
                     default.font = "Times New Roman", node.font.size = 20, edge.font.size = 8,
                     edge.line.style = 'SOLID',
                     edge.opacity = 175, edge.label.opacity = 255, border.opacity = 255, node.label.opacity = 255, node.fill.opacity = 255
                     ){

  if(!requireNamespace("RCy3", quietly = TRUE)){                                                                       # check RCy3 downloaded
    stop("In order to use this function, please download RCy3 as described in the vignette, the readme, and the function documentation.")
  }

  RCy3::cytoscapePing()                                                                                                 # ensure connection?

  # ACTUAL CODE AND DATA PROCESSING

  genes <- unique(c(cfn$Gene.1, cfn$Gene.2))                                                                            # get unique names of genes

  ptmnew <- ptmtable[, c("HCC4006_Erlotinib")]                                                                          # just take ptms and another col doesn't matter bc we replace
  rownames(ptmnew) <- sapply(rownames(ptmtable), function(x) strsplit(x, split = " ")[[1]][1])                          # get just the gene names in the PTM column
  rownames(ptmnew) <- c("score")                                                                                        # rename col

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
  cfn.nodes$node.type <- as.character(sapply(cfn.nodes$id, function(x) funckey$nodeType[which(funckey$Gene.Name == x)]))     # steal node type from funckey
  cfn.nodes$score <- as.numeric(sapply(cfn.nodes$id, function(x) sum(ptmnew$score[which(ptmnew$Gene == x)])))                # steal score from ptmnew

  cyscape <- RCy3::createNetworkFromDataFrames(cfn.nodes, cfn.edges, title = Network.title, collection = Network.collection) # create network (not sure if storing it does anything?)

  RCy3::copyVisualStyle("default", visual.style.name)   # create visual style



  # CUSTOMIZATION FROM HERE ON OUT

  RCy3::setNodeLabelMapping("id", style.name = visual.style.name)   # make the label names appear

  SetStandards(visual.style.name,                                   # we got standards in this joint
               background.color, edge.label.color, node.label.color, default.font, node.font.size, edge.font.size,
               edge.line.style,edge.opacity, edge.label.opacity, border.opacity, node.label.opacity, node.fill.opacity)

  NodeAppMap(visual.style.name)                                     # map the node colors CURRENTLY NOT SIZES

  EdgeAppMap(visual.style.name)                                     # map edge width CURRENTLY NOT COLORS

  NodeBorderMapping(visual.style.name)                              # map node border (color, thickness)

  RCy3::setVisualStyle(visual.style.name)                           # set vis style

}



# Graph CoCluster Correlation Network
#
# Creates a cytoscape graph of the cocluster correlation network. Ensure that you have the cytoscape app open and the RCy3 package downloaded and libraried.
#
# If error occurs, run cytoscapePing() to ensure connection to the Cytoscape Interface and try again.
#
# @param ptm.cccn A matrix showing strength of relationships between post-translational modifications using common clusters between the three distance metrics
# @param ptmtable A dataset for post-translational modifications. Formatted with numbered rows, and the first column containing PTM names. The rest of the column names should be drugs. Values are numeric values that represent how much the PTM has reacted to the drug.
# @param funckey A table graphing gene names to type of protein; defaults to internal database at PTMsToPathways::ex.funckey
# @param Network.title Desired title for the created Cytoscape Network; defaults to "cfn"
# @param Network.collection Desired name for the collection created on Cytoscape in which the network will reside; defaults to "PTMsToPathways"
# @param visual.style.name Desired name for the visual style created on Cytoscape; defaults to "PTMsToPathways.style"
#
# @param background.color Hex code of background color of graph; defaults to '#fcf3cf'
# @param edge.label.color Hex code of edge label color of graph; defaults to '#17202a'
# @param node.label.color Hex code of node label color of graph; defaults to '#145a32'
#
# @param default.font Font style of edge and node names; defaults to "Times New Roman"
# @param node.font.size Font size of the node name; defaults to 12
# @param edge.font.size Font size of the edge name; defaults to 8
#
# @param edge.line.style Type of edge style; defaults to "SOLID"; options include: "PARALLEL_LINES", "MARQUEE_EQUAL", "DOT", "EQUAL_DASH", "LONG_DASH", "CONTIGUOUS_ARROW", "MARQUEE_DASH", "DASH_DOT", "BACKWARD_SLASH", "FORWARD_SLASH", "VERTICAL_SLASH", "SOLID", "SEPARATE_ARROW", "MARQUEE_DASH_DOT", "ZIGZAG", "SINEWAVE"
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
#GraphPTMCccn <- function(ptm.cccn, ptmtable, funckey = PTMsToPathways::ex.funckey, Network.title = "cfn", Network.collection = "PTMsToPathways", visual.style.name = "PTMsToPathways.style",
#                     background.color = "#949494", edge.label.color = '#17202a', node.label.color = '#000000',
#                     default.font = "Times New Roman", node.font.size = 20, edge.font.size = 8,
#                     edge.line.style = 'SOLID',
#                     edge.opacity = 175, edge.label.opacity = 255, border.opacity = 255, node.label.opacity = 255, node.fill.opacity = 255
#                     ){
#
#  if(!requireNamespace("RCy3", quietly = TRUE)){                                                                       # check RCy3 downloaded
#    stop("In order to use this function, please download RCy3 as described in the vignette, the readme, and the function documentation.")
#    }
#
#  RCy3::cytoscapePing()                                                                                                 # ensure connection?
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
#}
