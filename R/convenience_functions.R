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
    both.seen <- apply(data,c(2,3),sum)
    X = both.seen
  }
  if(data_format == "GBI"){
    X = matrix(nrow = ncol(data), ncol = ncol(data))
    for(i in 1:nrow(X)){
      for(j in 1:nrow(X)){
        X[i,j] = sum(data[,i]==1 & data[,j]==1)
      }
    }
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
  denominator = matrix(nrow = dim(data)[2], ncol = dim(data)[2])
  if(!index %in% c("SRI", "HWI", "BII")) stop("Invalid Association Index")
  if(index == "BII" & data_format == "GBI") stop("BII not valid without sampling periods")
  if(data_format == "SP"){
    for(i in 1:nrow(denominator)){
      for(k in 1:nrow(denominator)){
        if(i > k){
          data.i = ifelse(rowSums(data[,i,]) > 0, 1, 0)
          data.k = ifelse(rowSums(data[,k,]) > 0, 1, 0)
          data.ik = data[,i,k]
          x = sum(data.ik)
          ya = sum(data.k == 0 & data.i == 1) #counts of only a
          yb = sum(data.i == 0 & data.k == 1) #counts of only b
          yab = sum(data.i == 1 & data.k == 1 & data.ik == 0) #both present, not seen together
          if(index == "HWI") denominator[i,k] = yab + 0.5 * (ya + yb) + x #HWI denominator
          if(index == "SRI") denominator[i,k] = yab + ya + yb + x #SRI denominator
          if(index == "BII") denominator[i,k] = yab + x #"both identified" association index
        }
      }
    }
  }
  if(data_format =="GBI"){
    for(i in 1:nrow(denominator)){
      for(j in 1:ncol(denominator)){
        x = sum(data[,i]==1 & data[,j]==1)
        ya = sum(data[,i]==1 & data[,j]==0)
        yb = sum(data[,i]==0 & data[,j]==1)
        if(index == "HWI") denominator[i,j] = 0.5*(ya + yb) + x
        if(index == "SRI") denominator[i,j] = ya+yb+x
      }
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
