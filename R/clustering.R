#' Hierarchical clustering of association indices using modularity
#'
#' An implementation of the network clustering method described by Lusseau (2007) and implemented in SOCPROG.
#'
#' @param network A symmetric square matrix of association indices
#' @param method The linkage method for hierarchical clustering
#'
#' @details This method works by transforming the association indices into distances (by subtracting them from 1), running hierarchical clustering as usual, and then calculating the modularity of the division at each split.
#' This method is currently only designed to be used with association-like indices (which have a maximum of 1).
#'
#' @return A named list containing the best modularity ("maximum_modularity"), the cophenetic clustering coefficient ("CCC"), the clustering tree ("tree"), a dataframe of cut heights and modularity ("cuts"), and merge information ("merge").
#' See \code{hclust} for details on manipulating and plotting trees.
#' @export

association_hclust <- function(network, method = "average"){
  m <- network
  graph = igraph::graph.adjacency(m, mode = "undirected", weighted = T)
  m.dist = stats::as.dist(1 - m)
  clustering = stats::hclust(m.dist, method = method)
  cp = stats::cophenetic(clustering)
  ccc = stats::cor(cp, m.dist)
  membership <- stats::cutree(clustering, k = 1:ncol(m))
  modularity <- apply(membership, 1, function(z){
    igraph::modularity(graph, z, weights = igraph::E(graph)$weight)
  })
  best_clusters = stats::cutree(clustering, k = which.max(modularity))
  height = 1 - cuts
  list(modularity = maximum_modularity, CCC = ccc, membership = best_clusters, tree = clustering, cuts = data.frame(Groups = 1:ncol(m), Mod = modularity), merge = clustering$merge)
}
