########## FIX THIS ONE IDK WHAT IS REALLY HAPPENING WHEN THE EXAMPLE IS CALLED ##########

#' Make GM Input Example
#' @keywords internal
ex.MakeDBInput <- function(ex.nodenames) {
  cat(ex.nodenames[[1]], sep = "\n")
}





#' Process GeneMANIA example
#' @keywords internal
ex.ProcessGMEdgefile <- function(epath, npath, ex.nodenames, ex.gm.network){
  gm.ex <- utils::head(ex.gm.network)
  print(gm.ex)
}





#' Find PPI Edges Example
#' @keywords internal
ex.FindPPIEdges <- function(cccn_matrix, db_filepaths = c(), gm.network = NA, ppi.network.name) {
  print(ppi.network.name)
}
