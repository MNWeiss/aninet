#' Generalised Linear Model Quadratic Assignment Procedure
#'
#' This function fits a generalised linear model of the given family
#' and assesses the significance of the estimated coefficients using a quadratic assignment procedure.
#' Multiple types of permutation are provided, with the default being the double semi-partialling method.
#'
#' @param formula A \code{glm} style formula describing the model to be fit. All elements of this formula should be square matrices of the same size.
#' @param family Either a \code{glm()} family or one of \code{"betar"} or \code{"negbin"}. See Details.
#' @param weights Square matrix of weights used for fitting, see Details.
#' @param offset Square matrix of offsets for fitting, see Details.
#' @param nperm Numeric, number of permutations to perform.
#' @param permutation Character specifying what form of permutation to carry out. One of \code{"DSP"}, \code{"Y"}, or \code{"X"}.
#'
#' @details The multiple regression quadratic assignment procedure (MRQAP) is the canonical method for dyadic regression in social networks. However, this method comes with two major assumptions: The residuals of the response are approximately normal (although the DSP method is fairly robust to violations of this assumption), and the social relationships of all dyads are measured with the same precision.
#' These assumptions are almost never met in animal social network analyses. Replacing the ordinary least-squares fit used in MRQAP with a GLM allows us to address both of these issues in a theoretically sound way, by specifying error families, offsets, and weights.
#'
#' The \code{weights} argument allows for specification of sampling weights per dyad. For binomial models (appropriate for association indices), these weights should be the dyadic denominator of the association index.
#'
#' The \code{offset} argument gives a matrix with a known coefficient of 1 in the model. This will primarily be useful for interaction rates. Here, we can use a count model without an upper bound (by setting \code{family} to be one of \code{"poisson"}, \code{"quasipoisson"}, or \code{"negbin"}), and include the logarithm of dyadic sampling effort as an offset. This means we're using sampling effort as an exposure term, and therefore modelling interaction rates rather than just counts.
#'
#' In most cases, the model is fit using \code{glm()}. However, if \code{family = "betar"}, the \code{betareg::betareg()} function is used, and if \code{family = "negbin"}, the \code{MASS::glm.nb()} function is used.
#'
#' Beta models will be most useful for association index type data without an integer numerator/denominator (such as measurements of portion time together from biologgers/video). Dyadic sampling effort can still be specified as the weights, but take care.
#' In \code{betareg()}, weights are treated as sampling, rather than proportional, weights. This means that the function assumes your true sample size is \code{sum(weights)}. While this won't effect your estimate or significance (because we use permutations), it will give
#' pretty strange results for the standard errors. A solution is to transform your weights such that \code{sum(weights) = length(weights)}.
#'
#' The function allows multiple types of permutation. The \code{"DSP"} method is the most robust for testing multiple predictors, and is based on the method proposed by Dekker et al. (2007). This explicitly tests the effect of each covariate, controlling for the effect of others and the relationship between variables.
#' The \code{"Y"} permutation method permutes the response matrix, and tests the null hypothesis that the response is unrelated to any of the predictors.
#' The \code{"X"} method permutes each predictor matrix, testing the null hypothesis of no relationship, but importantly does not control for any covariance among predictors.
#' Note that if your formula contains only a single predictor, the function will default to \code{permutation = "Y"}.
#'
#' @return An object of class \code{glmqap}, containing a summary of the fitted model. The coefficients table now contains p-values from the randomizations. In addition to all information normally in the respective model summary objects, the object contains an element called \code{permuted_z}, which is a matrix containing the permuted value of each pivotal statistic for all permutations.
#'
#' @export
glmqap <- function(formula, family = "gaussian", weights=NULL, offset=NULL, nperm=1000, permutation = "DSP"){

  fm <- formula

  formula <- stats::model.frame(formula, na.action = NULL)
  i <- attr(attr(formula, "terms"), "response")
  response <- as.matrix(formula[i])
  x_names <- attr(attr(formula, "terms"), "term.labels")
  predictors <- list()

  for (i in 1:length(x_names)) {
    predictors[[i]] <- as.matrix(formula[[x_names[i]]])
  }

  if(!all(unlist(lapply(predictors, function(z) nrow(z) == ncol(z) )))){
    stop("Predictors must be square matrices")
  }
  if(nrow(response) != ncol(response)){
    stop("Response must be a square matrix")
  }

  if(!all(unlist(lapply(predictors, function(z) nrow(z) == nrow(response) )))){
    stop("Predictors and response must be of the same dimensions")
  }

  if(family == "binomial" & any(!is.integer(response)) & is.null(weights)){
    stop("Provide denominators as weights for binomial model")
  }

  if(!is.null(weights)){
    if(!is.matrix(weights) | nrow(weights) != ncol(weights) | nrow(weights) != nrow(response)){
      stop("Weights must be a square matrix with same dimensions as response")
    }
  }
  if(!is.null(offset)){
    if(!is.matrix(offset) | nrow(offset) != ncol(offset) | nrow(offset) != nrow(response)){
      stop("Offsets must be a square matrix with same dimensions as response")
    }
  }

  names(predictors) <- x_names

  if(length(x_names) == 1) permutation <- "Y"

  y <- response[lower.tri(response)] #vectorized response
  x <- do.call(cbind,lapply(predictors, function(z) z[lower.tri(z)])) #matrix of predictors

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

  message("Fitting Model")

  if(family != "betar" & family != "negbin"){ #if a built-in GLM family is specified
    mod.orig <- stats::glm(y ~ x, weights = w, offset = o, family = family) #fit the GLM
    z.val <- summary(mod.orig)$coefficients[-1,3] #pull out predictor z values
    coef <- mod.orig$coefficients
    se <- summary(mod.orig)$coefficients[,2]
  }
  if(family == "betar"){ #if beta regression is specified
    mod.orig <- betareg::betareg(y ~ x, weights = w) #fit model
    z.val <- summary(mod.orig)$coefficients$mean[-1,3] #pull out z values
    coef <- mod.orig$coefficients$mean
    se <- summary(mod.orig)$coefficients$mean[,2]
  }
  if(family == "negbin"){ #same for negative binomial regression
    mod.orig <- MASS::glm.nb(y ~ x + offset(o), weights = w)
    z.val <- summary(mod.orig)$coefficients[-1,3]
    coef <- mod.orig$coefficients
    se <- summary(mod.orig)$coefficients[,2]
  }

  summ <- summary(mod.orig)
  summ$call[[2]] <- fm
  summ$call <- summ$call[1:2]

  z.perm <- matrix(nrow = nperm, ncol = length(z.val)) #matrix to hold permuted z values

  message("Performing Permutations")

  if(permutation == "DSP"){

    pb <- utils::txtProgressBar(min = 1, max = nperm*length(predictors), style = 3) #create a progress bar
    counter <- 1 #and a counter

    for(i in 1:length(predictors)){ #for each predictor

      x <- predictors[[i]] #pull out the predictor of interest
      z <- predictors[-i] #other predictors

      x <- x[lower.tri(x)] #vectorize
      z <- do.call(cbind,lapply(z,function(p)p[lower.tri(p)])) #matrix

      eps <- stats::residuals(stats::lm(x ~ z)) #residuals of x given z

      eps.mat <- response #turn it into a matrix
      eps.mat[lower.tri(eps.mat)] <- eps
      eps.mat[upper.tri(eps.mat)] <- t(eps.mat)[upper.tri(eps.mat)] #symmetrize

      for(j in 1:nperm){

        utils::setTxtProgressBar(pb,counter) #update progress
        counter <- counter+1 #update counter

        samp.order <- sample(ncol(response)) #permuted order of rows and columns
        eps.perm <- eps.mat[samp.order,samp.order] #permuted residuals
        eps.perm <- eps.perm[lower.tri(eps.perm)] #vectorized

        #fit permuted model and pull out z statistic

        if(family != "betar" & family != "negbin"){
          z.perm[j,i] <- summary(stats::glm(y ~ eps.perm + z, family = family, offset = o, weights = w))$coefficients[2,3]
        }
        if(family == "betar"){
          z.perm[j,i] <- summary(betareg::betareg(y ~ eps.perm + z, weights = w))$coefficients$mean[2,3]
        }
        if(family == "negbin"){
          z.perm[j,i] <- summary(MASS::glm.nb(y ~ eps.perm + z + offset(o), weights = w))$coefficients[2,3]
        }

      }

    }

  }

  if(permutation == "X"){

    pb <- utils::txtProgressBar(min = 1, max = nperm*length(predictors), style = 3) #create a progress bar
    counter <- 1 #and a counter

    for(i in 1:length(predictors)){ #for each predictor

      x <- predictors[[i]] #pull out the predictor of interest
      z <- predictors[-i] #other predictors

      z <- do.call(cbind,lapply(z,function(p)p[lower.tri(p)])) #matrix

      for(j in 1:nperm){

        utils::setTxtProgressBar(pb,counter) #update progress
        counter <- counter+1 #update counter

        samp.order <- sample(ncol(response)) #permuted order of rows and columns
        xp <- x[samp.order,samp.order]
        xp <- xp[lower.tri(xp)]

        #fit permuted model and pull out z statistic

        if(family != "betar" & family != "negbin"){
          z.perm[j,i] <- summary(stats::glm(y ~ xp + z, family = family, offset = o, weights = w))$coefficients[2,3]
        }
        if(family == "betar"){
          z.perm[j,i] <- summary(betareg::betareg(y ~ xp + z, weights = w))$coefficients$mean[2,3]
        }
        if(family == "negbin"){
          z.perm[j,i] <- summary(MASS::glm.nb(y ~ xp + z + offset(o), weights = w))$coefficients[2,3]
        }

      }

    }

  }

  if(permutation == "Y"){
    pb <- utils::txtProgressBar(min = 1, max = nperm, style = 3) #create a progress bar
    counter <- 1 #and a counter
    for(i in 1:nperm){
      utils::setTxtProgressBar(pb,counter) #update progress
      counter <- counter+1 #update counter
      samp_order <- sample(nrow(response))
      r.p <- response[samp_order,samp_order]
      yp <- r.p[lower.tri(r.p)]

      if(!is.null(offset)){
        offset.p <- offset[samp_order,samp_order]
        op <- offset.p[lower.tri(offset.p)]
      }else{
        op <- NULL
      }
      if(!is.null(weights)){
        weights.p <- weights[samp_order,samp_order]
        wp <- weights.p[lower.tri(weights.p)]
      }else{
        wp <- NULL
      }

      if(family != "betar" & family != "negbin"){
        z.perm[i,] <- summary(stats::glm(yp ~ x, family = family, offset = op, weights = wp))$coefficients[-1,3]
      }
      if(family == "betar"){
        z.perm[i,] <- summary(betareg::betareg(yp ~ x, weights = wp))$coefficients$mean[-1,3]
      }
      if(family == "negbin"){
        z.perm[i,] <- summary(MASS::glm.nb(yp ~ x + offset(op), weights = wp))$coefficients[-1,3]
      }

    }

  }

  pval <- sapply(1:length(predictors), function(p){
    min(
      c(sum(z.perm[,p] >= z.val[p])+1,
        sum(z.perm[,p] <= z.val[p])+1
      ))/(nperm+1)
  })*2

  res <- list(
    call = fm,
    pred_names = x_names,
    family = family,
    permutation = permutation,
    nperm = nperm,
    coefficients = coef,
    stderr = se,
    z = c(NA,z.val),
    p = c(NA,pval),
    permuted_z = z.perm,
    aic = AIC(mod.orig),
    bic = BIC(mod.orig),
    loglik = as.numeric(logLik(mod.orig)),
    weights = deparse(substitute(weights)),
    offset = deparse(substitute(offset))
  )

  if(family == "betar"){
    res$precision <- mod.orig$coefficients$precision
  }
  if(family == "negbin"){
    res$theta <- mod.orig$theta
  }

  class(res) <- "glmqap"

  return(res)

}

#' Print function for glmqap objects
#'
#' Prints summary information about GLMQAP models
#'
#' @param x A glmqap object
#' @details Prints formatted results from GLMQAP fits, including formula, family, permutations, weights and offsets, coefficients, standard errors, test statistics, and p-values
#' @return Printed results
#' @method print glmqap
#' @export
print.glmqap <- function(x){

  perm_method <- ifelse(x$permutation == "X", "X permutation", ifelse(x$permutation == "Y", "Y permutation", "Doube-Semi-Partialling"))
  cat(paste("GLMQAP with ", perm_method, "\n\n", sep = ""))
  cat(paste("Formula: ", Reduce(paste, deparse(x$call)), "\n"))
  cat(paste("Family: ", x$family), "\n")
  cat(paste("Weights:", x$weights), "\n")
  cat(paste("Offset:", x$offset), "\n")
  cat(paste("Permutations: ", x$nperm), "\n\n")
  cat("Coefficients:\n")
  results_mat <- cbind(x$coefficients, x$stderr, x$z, x$p)
  colnames(results_mat) <- c("Estimate", "Std. Error", "Z", "P(two-tailed)")
  row.names(results_mat) <- c("Intercept",x$pred_names)
  print.table(results_mat)
  cat("\n")
  cat("log-likelihood:", x$loglik, "\t")
  cat("AIC:", x$aic, "\t")
  cat("BIC:", x$bic, "\t")
  cat("\n")

}
