#' Generalised Affiliation Indices
#'
#' Caluclates a social network based on generalised affiliation indices (GAI). When non-social factors influence the observed interactions or associations, using GAIs can help to remove the effect of these factors, revealing underlying social affinities. This is done by fitting a generalised linear model, with the observed network as a response and non-social confounds as predictors, and using the residuals as a measure of affiliation.
#'
#' @param formula A glm style formula giving the response and predictor matrices.
#' @param family Error family for fitting, either a \code{glm} family or one of \code{"betar"} or \code{"negbin"}.
#' @param weights Weights used for fitting, see Details.
#' @param offset Offset used in fitting, see Details.
#' @param type Type of residuals to calculate
#'
#' @details For association data, \code{family = "binomial"} is the most appropriate. In this case, the observed association indices should be the response matrix, and the dyadic denominators should be passed to the \code{weights} argument.
#' For interaction rates, one of \code{"poisson"}, \code{"quasipoisson"}, or \code{"negbin"} should be passed to \code{family}, with log(sampling effort) passed as an offset.
#' For completeness, beta regression is also included, using \code{betareg::betareg()}.
#' Generally, its a good idea to test whether the predictors are important prior to calculating GAIs. This can be done using GLMQAP.
#'
#' @return A square matrix containing the generalised affiliation indices.
#' @export
gai <- function(formula,family="gaussian",weights=NULL,offset=NULL,type=c("deviance","pearson","response")){

  if(length(type) > 1){
    type <- type[1]
  }

  fm <- formula

  formula <- stats::model.frame(formula, na.action = NULL)
  i <- attr(attr(formula, "terms"), "response")
  response <- as.matrix(formula[i])
  x_names <- attr(attr(formula, "terms"), "term.labels")
  predictors <- list()

  for (i in 1:length(x_names)) {
    predictors[[i]] <- as.matrix(formula[[x_names[i]]])
  }

  names(predictors) <- x_names

  y <- response[lower.tri(response)] #vectorized response
  x <- do.call(cbind,lapply(predictors,function(z)z[lower.tri(z)])) #matrix of predictors

  if(is.null(weights)){ #if no weights are specified
    w <- NULL #w is null
  }else{ #otherwise
    w <- weights[lower.tri(weights)] #pull out weights as a vector
  }

  if(is.null(offset)){ #same as above for offsets
    o <- NULL
  }else{
    o <- offset[lower.tri(offset)]
  }

  if(family != "betar" & family != "negbin"){ #if a built-in GLM family is specified
    mod.orig <- stats::glm(y ~ x, weights = w, offset = o, family = family) #fit the GLM
  }
  if(family == "betar"){ #if beta regression is specified
    mod.orig <- betareg::betareg(y ~ x, weights = w) #fit model
  }
  if(family == "negbin"){ #same for negative binomial regression
    mod.orig <- MASS::glm.nb(y ~ x + offset(o), weights = w)
  }

  resid <- stats::residuals(mod.orig,type=type)

  gai.mat <- response
  gai.mat[lower.tri(gai.mat)] <- resid
  gai.mat[upper.tri(gai.mat)] <- t(gai.mat)[upper.tri(gai.mat)]
  diag(gai.mat) <- 0

  colnames(gai.mat) <- row.names(gai.mat)

  return(gai.mat)

}

#' Attribute similarity matrices
#'
#' Generate square matrices of similarity between individuals from individual attribute data
#'
#' @param attr Vector, individual attributes
#' @param type Character, type of similarity to calculate
#'
#' @details The available types of similarity are \code{"discrete"} for 1-0 similarity, \code{"absdiff"} for negative absolute difference, and \code{"sqdiff"} for negative squared difference. If the attribute data is not numeric, the function will default to \code{"discrete"}.
#' @return A square matrix of similarity values
#' @export
attribute_similarity <- function(attr,type="discrete"){
  if(class(attr) != "numeric" & class(attr) != "integer") type = "discrete"
  if(type == "discrete"){
    res <- sapply(attr,function(z)ifelse(attr==z,1,0))
  }
  if(type == "absdiff"){
    res <- sapply(attr,function(z) - abs(attr-z))
  }
  if(type == "sqdiff"){
    res <- sapply(attr,function(z) - (attr-z)^2)
  }
  diag(res) <- 0
  res
}

#' Joint gregariousness
#'
#' This function calculates the "joint gregariousness" for each dyad, as given by Whitehead & James (2015)
#'
#' @param network A square matrix, defining the network from which to calculate joint gregariousness
#' @param log.jg Logical, indicating whether the log-transformed scores are to be returned.
#'
#' @details Joint gregariousness is defined by Whitehead & James (2015) as:
#' \deqn{joint gregariousness[i,j] = log(greg[i] * greg[j])}, where \eqn{greg[i]} is \deqn{\sum_{k \neq j} network[i,k]}
#' This definition removes each dyads edge weight from the calculation of their joint gregariousness.
#' The function allows the user to get the non log-transformed value. This is useful if there are individuals with only as single non-zero edge weight.
#'
#' @return A square matrix of dyadic joint gregariousness scores.
#' @export
joint_gregariousness <- function(network,log.jg=T){

  j.g = network

  for(i in 1:nrow(j.g)){
    for(k in 1:nrow(j.g)){
      greg.i = sum(network[i, -k])
      greg.k = sum(network[k, -i])
      if(log.jg){
        j.g[i,k] = log(greg.i * greg.k)
      }else{
        j.g[i,k] = greg.i * greg.k
      }
    }
  }

  diag(j.g) = 0

  j.g

}
