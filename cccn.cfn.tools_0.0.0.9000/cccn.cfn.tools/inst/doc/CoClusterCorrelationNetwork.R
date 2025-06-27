## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "plots/",
  fig.width = 8,
  fig.height = 6,
  out.width = "100%",
  results = "hold"
)

## ----eval = FALSE-------------------------------------------------------------
# MakeClusterList(ptmtable, correlation.matrix.name = "ptm.correlation.matrix", list.name = "tsne.matrices", toolong = 3.5)

## ----echo = FALSE, eval = TRUE------------------------------------------------
load('vig_figs/eu-clusters.rda')
head(eu_ptms_list, 1)

## ----eval = FALSE-------------------------------------------------------------
# MakeCorrelationNetwork(tsne.matrices, ptm.correlation.matrix, keeplength = 2, common.clusters.name = "common.clusters", cccn.name = "cccn_matrix")

## ----eval = FALSE-------------------------------------------------------------
# MakeDBInput(cccn_matrix, file.path.name = "db_nodes.txt")

## ----eval = FALSE-------------------------------------------------------------
# ProcessGMEdgefile(gm.edgefile.path, gm.nodetable.path, nodenames, gm.network.name = "gm.network")

## ----eval = FALSE-------------------------------------------------------------
# FindPPIEdges(cccn_matrix, db_file_paths = c(), ppi.network.name = "ppi.network")

## ----eval = FALSE-------------------------------------------------------------
# ClusterFilteredNetwork(cccn_matrix, ppi.network, cfn.name = "cfn")

