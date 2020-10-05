#' Calculate composite sociality indices from a set of social networks
#'
#' Collapses multiple types of interaction or association into a single network of composite indices of social relationships
#'
#' @param networks A list of networks from which to calculate composite indices. Must be of the same dimensions, and in the same order.
#'
#' @details Dyadic composite sociality indices are calculated as in Sapolsky et al. (1991).
#' The dyadic CSI is calculated as:
#'
#' \deqn{CSI_{ij} = frac{sum_{d = 1}^{D} frac{f_{dij}}{\bar{f_d}} } {D}}
#'
#' Before calculating these indices, it is important to consider whether the included networks positively covary. For example, grooming and proximity may be useful for generating a composite index, but it likely does not make sense to combine grooming and aggression.
#'
#' @return A square matrix of dyadic compisite sociality indices
#'
#' @export

dyadic_csi <- function(networks){

  nrows <- unlist(lapply(networks, nrow))
  ncols <- unlist(lapply(networks, ncol))
  classes <- unlist(lapply(networks, class))
  if(length(unique(nrows)) > 1) stop("Networks not of same dimensions!")
  if(any(nrows != ncols)) stop("Networks must be square matrices")
  if(any(classes != "matrix")) stop("Networks must me matrices")

  n <- unique(nrows)
  d <- length(networks)

  mean_rate <- unlist(lapply(networks,function(z) mean(z[lower.tri(z)|upper.tri(z)]) ))
  std_rate <- lapply(networks, function(z){
    mu <- mean(z[lower.tri(z)|upper.tri(z)])
    r <- z/mu
    r
  })
  net_array <- array(dim = c(d,n,n))
  for(i in 1:d){
    net_array[i,,] <- std_rate[[i]]
  }

  csi <- apply(net_array,c(2,3),mean)

  return(csi)

}
