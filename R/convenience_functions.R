#' Get numerator of association indices
#'
#' For many functions, the numerator of association indices is needed, rather than the final network. This function allows to user to easily extract these.
#'
#' @param data The association data
#' @param data_format Character, indicating what form the data is in.
#' @param return Character, indicating whether a vector or square matrix should be returned
#'
#' @details The association index does not need to be specified, as they almost all use the same numerator. Data format can be either \code{"SP"} or \code{"GBI"}, the two main data formats used by \code{"asnipe"}.
#' @return Either a vector of numerators (if \code{return = "vector"}) or a square matrix of numerators (if \code{return = "matrix"}).
#' @export
get_numerator <- function(data, data_format = "SP", return = "vector"){
  if(data_format == "SP"){
    X <- apply(data,c(2,3),sum)
  }
  if(data_format == "GBI"){
    X <- t(data) %*% data
    diag(X) <- 0
  }
  if(return == "vector" | is.null(return)){
    return(X[lower.tri(X)])
  }
  if(return == "matrix"){
    diag(X) = 0
    colnames(X) = colnames(data)
    row.names(X) = colnames(data)
    return(X)
  }
}

#' Get denominator of association indices
#'
#' For many functions, the denominator of association indices is needed, rather than the final network. This function allows to user to easily extract these.
#'
#' @param data The association data.
#' @param index The association index to use, one of \code{"SRI"}, \code{"HWI"}, and \code{"BII"}.
#' @param data_format Character, indicating what form the data is in.
#' @param return Character, indicating whether a vector or square matrix should be returned
#'
#' @details Data format can be either \code{"SP"} or \code{"GBI"}, the two main data formats used by \code{"asnipe"}. Passing \code{index = "SRI"} returns the denominator of the simple ratio index, \code{index = "HWI"} returns half-weight index, while \code{index = "BII"} uses a "both-identified" index.
#' @return Either a vector of denominators (if \code{return = "vector"}) or a square matrix of denominators (if \code{return = "matrix"}).
#' @export
get_denominator <- function(data, index = "SRI", data_format = "SP", return = "vector"){
  if(!index %in% c("SRI", "HWI", "BII")) stop("Invalid Association Index")
  if(data_format == "SP"){
    occur <- ifelse(apply(data,c(1,2),sum) > 0, 1, 0)
    if(index == "BII") denominator <- t(occur) %*% occur
    if(index == "SRI") denominator <- nrow(occur) - t(1-occur) %*% (1-occur)
    if(index == "HWI"){
      ya <- t(occur) %*% (1-occur)
      yb <- t(ya)
      yab <- t(occur) %*% (occur)
      denominator <- 0.5*(ya+yb) + yab
    }
  }
  if(data_format =="GBI"){
    if(index == "BII") denominator <- t(data) %*% data
    if(index == "SRI") denominator <- nrow(data) - t(1-data) %*% (1- data)
    if(index == "HWI"){
      ya <- t(data) %*% (1-data)
      yb <- t(ya)
      yab <- t(data) %*% data
      denominator <- 0.5*(ya+yb) + yab
    }
  }
  if(return == "vector" | is.null(return)){
    return(denominator[lower.tri(denominator)])
  }
  if(return == "matrix"){
    D = denominator
    D[upper.tri(D)] = t(D)[upper.tri(D)]
    colnames(D) = colnames(data)
    row.names(D) = colnames(data)
    diag(D) = 0
    return(D)
  }
}
