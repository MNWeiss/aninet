#' Latent social space models for interaction and association data
#'
#' Fits a D dimensional latent space model of the given family to a matrix of associations or interactions
#'
#' @param formula A glm style formula, where the response and predictors are all square matrices of the same dimension. The response matrix should be the counts of interactions or associations.
#' @param family Character, one of either "poisson" or "binomial". See details.
#' @param dimensions Integer, the number of latent dimensions to model.
#' @param ind.RE Logigal, indicating whether to include an individual sociality random effect.
#' @param effort A square matrix indicating dyadic sampling effort. See details.
#' @param beta.prior Numeric vector, prior for fixed effects. Should be the mean and precision of a normal distribution.
#' @param vcv.prior Numeric matrix, prior for the variance-covariance matrix of latent positions. Should be a matrix with row and column number equal to the dimensions, representing a parameterization of the Wishart distribution.
#' @param re.prior Numeric vector, prior for individual random effect. Should be the parameters of a Gamma distribution.
#' @param z0 Optional numeric matrix, initial positions of nodes in the latent space. If not provided, initial values are generated using mutlidimensional scaling.
#' @param ... Further arguments to be passed to runjags.
#'
#' @details Social networks often exhibit transitivity, where a connection between B and C is likely to be stronger if B and C have strong connections to a third individual A. This tendency is often referred to as "triadic closure."
#' One way to account for this in regression settings is to view nodes as being placed within a D dimensional latent space, with edges partially determined by the euclidean distances between nodes. Due to the triangle inequality, the distance between B and C (D(B,C)) will always be less than or equal to the sum of the distances D(A,C) and D(A,B). Thus, the model induces some degree of transitivity.
#' This function fits one of these latent space models to a matrix of associations or interactions. For interactions, a Poisson model should be fit, with effort indicating the sampling time per dyad. For associations, a binomial model should be fit, with effort being indicated by the denominator of the association index.
#' In both cases, the response matrix should be a matrix of integers, indicating the number of dyadic interactions or associations.
#' This function fits the model using Gibbs sampling via JAGS and runjags, and allows the user to define priors for fixed effects, the covariance matrix of the latent positions, and individual random effects (if included).
#'
#'
#'
#' @return A named list containing the following slots
#'
#' @slot summary The summary statistics for the fixed effects and variance of the random effects (if included).
#' @slot distances The sampled distances between nodes in social space, stored as a S x N x N array, where S is the number of MCMC samples.
#' @slot z The sampled positions of nodes in social space, stored as a S x N x D array, where S is the number of MCMC samples and D is the number of latent dimensions.
#' @slot jags_model The full runjags object, which can be used to assess convergence through trace plots and diagnostics, and to get the DIC value of the model.
#'
#' @export
latent_space <- function(formula, family = "poisson", dimensions = 2, ind.RE = T, effort, beta.prior = c(0,1e-4), vcv.prior = NULL, re.prior = c(0.1,0.1), z0 = NULL, ...){

  if(!is.integer(dimensions)) dimensions <- round(dimensions)
  if(dimensions == 0) stop("Latent space model must have at least 1 latent dimension")

  if(is.null(vcv.prior)){
    vcv.prior <- diag(1,dimensions)
  }

  if(class(vcv.prior)[1] != "matrix" | nrow(vcv.prior) != dimensions){
    stop("Prior for covariance matrix must have row and column number equal to latent space dimensions")
  }
  if(length(beta.prior) != 2 | !is.numeric(beta.prior)){
    stop("Prior for fixed effects should be the mean and precision of a normal distribution")
  }
  if(length(re.prior) != 2 | !is.numeric(re.prior)){
    stop("Prior for random effects should be the parameters of a Gamma distribution")
  }

  cat("Setting up model structure...")
  cat("\n")

  mf <- stats::model.frame(formula, na.action = NULL)
  i <- attr(attr(mf, "terms"), "response")
  response <- as.matrix(mf[i])

  if(!is.null(z0)){
    if(nrow(z0) != nrow(response) | ncol(z0) != dimensions){
      stop("Initial Z values of wrong dimensions")
    }
  }

  if(is.null(z0)){
    cat("No initial configuration provided, generating initial points using MDS...")
    cat("\n")
    rate <- response/effort
    diag(rate) <- 0
    rate.dist <- 1 - (rate/max(rate, na.rm=T)) + 1e-6
    diag(rate.dist) <- 0
    rate.dist <- as.dist(rate.dist)
    z0 <- MASS::isoMDS(rate.dist,k=dimensions)$points
  }

  x_names <- attr(attr(mf, "terms"), "term.labels")

  p <- length(x_names)
  n <- ncol(response)

  if(p > 0){

    predictors <- list()

    for (i in 1:length(x_names)) {
      predictors[[i]] <- as.matrix(mf[[x_names[i]]])
    }

    if(any(unlist(lapply(predictors,function(z)!is.matrix(z)|!isSymmetric(z)|nrow(z)!=nrow(response))))){
      stop("Predictors must be symmetric square matrices of same dimension as the response")
    }

  }

  pred_array <- array(dim = c((p+1),n,n))
  for(i in 1:(p+1)){
    if(i == 1){
      pred_array[i,,] <- 1
    }else{
      pred_array[i,,] <- predictors[[(i-1)]]
    }
  }

  if(ind.RE){
    data <- list(
      G = dimensions,
      N = n,
      effort = effort,
      P = p+1,
      zRmat = vcv.prior,
      y = response,
      x = pred_array,
      beta.prior = beta.prior,
      re.prior = re.prior
    )
  }else{
    data <- list(
      G = dimensions,
      N = n,
      effort = effort,
      P = p+1,
      zRmat = vcv.prior,
      y = response,
      x = pred_array,
      beta.prior = beta.prior
    )
  }

  if(family == "poisson"){

    if(ind.RE){

      filestring <- "model{

        for(i in 1:(N-1)){
          for(j in (i+1):N){

            y[i,j] ~ dpois(lambda[i,j])

            lambda[j,i] <- lambda[i,j]
            log(lambda[i,j]) <- sum(beta[1:P]*x[1:P,i,j]) - d[i,j] + re[i] + re[j] + log(effort[i,j])


            d[j,i] <- d[i,j] #symmetrize
            d[i,j] <- sqrt(sum((z[i,1:G]-z[j,1:G])^2)) #get the distances from estimates of z

          }
        }

      for(i in 1:P){
        beta[i] ~ dnorm(beta.prior[1],beta.prior[2])
      }

      for(i in 1:N){

        re[i] ~ dnorm(0,tau)

        d[i,i] <- 0
        lambda[i,i] <- 0

        z[i,1:G] ~ dmnorm(zMu, zInvCovMat[1:G,1:G])

      }

      tau <- pow(sigma,-2)
      sigma ~ dgamma(re.prior[1],re.prior[2]) #prior for individual variance

      for(i in 1:G){

        zMu[i] <- 0 #we'll center the latent space at the origin

      }

      zInvCovMat ~ dwish(zRmat[1:G,1:G], G) #prior for the spatial variance-covariance matrix

    }"

    params = c("beta","z","sigma","dic")

    }else{

      filestring <- "model{

        for(i in 1:(N-1)){
          for(j in (i+1):N){

            y[i,j] ~ dpois(lambda[i,j])

            lambda[j,i] <- lambda[i,j]
            log(lambda[i,j]) <- sum(beta[1:P]*x[1:P,i,j]) - d[i,j] + log(effort[i,j])


            d[j,i] <- d[i,j] #symmetrize
            d[i,j] <- sqrt(sum((z[i,1:G]-z[j,1:G])^2)) #get the distances from estimates of z

          }
        }

        for(i in 1:P){
          beta[i] ~ dnorm(beta.prior[1],beta.prior[2])
        }

        for(i in 1:N){

          d[i,i] <- 0
          lambda[i,i] <- 0

          z[i,1:G] ~ dmnorm(zMu, zInvCovMat[1:G,1:G])

        }

        for(i in 1:G){

          zMu[i] <- 0 #we'll center the latent space at the origin

        }

        zInvCovMat ~ dwish(zRmat[1:G,1:G], G) #prior for the spatial variance-covariance matrix

      }"

      params <- c("beta","z","dic")

    }

  }

  if(family == "binomial"){

    if(ind.RE){

      filestring <- "model{

        for(i in 1:(N-1)){
          for(j in (i+1):N){

            y[i,j] ~ dbinom(mu[i,j],effort[i,j])

            mu[j,i] <- mu[i,j]
            logit(mu[i,j]) <- sum(beta[1:P]*x[1:P,i,j]) - d[i,j] + re[i] + re[j]


            d[j,i] <- d[i,j] #symmetrize
            d[i,j] <- sqrt(sum((z[i,1:G]-z[j,1:G])^2)) #get the distances from estimates of z

          }
        }

      for(i in 1:P){
        beta[i] ~ dnorm(beta.prior[1],beta.prior[2])
      }

      for(i in 1:N){

        re[i] ~ dnorm(0,tau)

        d[i,i] <- 0
        lambda[i,i] <- 0

        z[i,1:G] ~ dmnorm(zMu, zInvCovMat[1:G,1:G])

      }

      tau <- pow(sigma,-2)
      sigma ~ dgamma(re.prior[1],re.prior[2]) #prior for individual variance

      for(i in 1:G){

        zMu[i] <- 0 #we'll center the latent space at the origin

      }

      zInvCovMat ~ dwish(zRmat[1:G,1:G], G) #prior for the spatial variance-covariance matrix

    }"

      params = c("beta","z","sigma","dic")

    }else{

      filestring <- "model{

        for(i in 1:(N-1)){
          for(j in (i+1):N){

            y[i,j] ~ dbinom(mu[i,j],effort[i,j])

            mu[j,i] <- mu[i,j]
            logit(mu[i,j]) <- sum(beta[1:P]*x[1:P,i,j]) - d[i,j]


            d[j,i] <- d[i,j] #symmetrize
            d[i,j] <- sqrt(sum((z[i,1:G]-z[j,1:G])^2)) #get the distances from estimates of z

          }
        }

        for(i in 1:P){
          beta[i] ~ dnorm(beta.prior[1],beta.prior[2])
        }

        for(i in 1:N){

          d[i,i] <- 0
          lambda[i,i] <- 0

          z[i,1:G] ~ dmnorm(zMu, zInvCovMat[1:G,1:G])

        }

        for(i in 1:G){

          zMu[i] <- 0 #we'll center the latent space at the origin

        }

        zInvCovMat ~ dwish(zRmat[1:G,1:G], G) #prior for the spatial variance-covariance matrix

      }"

      params <- c("beta","z","dic")

    }

  }

  if(family == "beta"){

    if(ind.RE){

      filestring <- "model{

        for(i in 1:(N-1)){
          for(j in (i+1):N){

            y[i,j] ~ dbinom(mu[i,j],effort[i,j])

            mu[j,i] <- mu[i,j]
            logit(mu[i,j]) <- sum(beta[1:P]*x[1:P,i,j]) - d[i,j] + re[i] + re[j]


            d[j,i] <- d[i,j] #symmetrize
            d[i,j] <- sqrt(sum((z[i,1:G]-z[j,1:G])^2)) #get the distances from estimates of z

          }
        }

      for(i in 1:P){
        beta[i] ~ dnorm(beta.prior[1],beta.prior[2])
      }

      for(i in 1:N){

        re[i] ~ dnorm(0,tau)

        d[i,i] <- 0
        lambda[i,i] <- 0

        z[i,1:G] ~ dmnorm(zMu, zInvCovMat[1:G,1:G])

      }

      tau <- pow(sigma,-2)
      sigma ~ dgamma(re.prior[1],re.prior[2]) #prior for individual variance

      for(i in 1:G){

        zMu[i] <- 0 #we'll center the latent space at the origin

      }

      zInvCovMat ~ dwish(zRmat[1:G,1:G], G) #prior for the spatial variance-covariance matrix

    }"

      params = c("beta","z","sigma","dic")

    }else{

      filestring <- "model{

        for(i in 1:(N-1)){
          for(j in (i+1):N){

            y[i,j] ~ dbinom(mu[i,j],effort[i,j])

            mu[j,i] <- mu[i,j]
            logit(mu[i,j]) <- sum(beta[1:P]*x[1:P,i,j]) - d[i,j]


            d[j,i] <- d[i,j] #symmetrize
            d[i,j] <- sqrt(sum((z[i,1:G]-z[j,1:G])^2)) #get the distances from estimates of z

          }
        }

        for(i in 1:P){
          beta[i] ~ dnorm(beta.prior[1],beta.prior[2])
        }

        for(i in 1:N){

          d[i,i] <- 0
          lambda[i,i] <- 0

          z[i,1:G] ~ dmnorm(zMu, zInvCovMat[1:G,1:G])

        }

        for(i in 1:G){

          zMu[i] <- 0 #we'll center the latent space at the origin

        }

        zInvCovMat ~ dwish(zRmat[1:G,1:G], G) #prior for the spatial variance-covariance matrix

      }"

      params <- c("beta","z","dic")

    }

  }

  fit <- runjags::run.jags(model = filestring, data = data, monitor = params, ...)

  res <- do.call(rbind,fit$mcmc)

  cat("Performing transformations on latent positions...")
  cat("\n")

  z <- array(dim = c(nrow(res), n, dimensions))

  for(i in 1:nrow(res)){
    zi <- matrix(res[i,substr(colnames(res),1,1) == "z"], nrow = n, ncol = dimensions, byrow = F)
    z[i,,] <- MCMCpack::procrustes(X = zi, Xstar = z0, translation = F, dilation = F)$X.new
  }

  d <- array(dim = c(nrow(res),n,n))
  for(i in 1:nrow(res)) d[i,,] <- as.matrix(dist(z[i,,]))

  if(ind.RE){
    summ <- summary(fit, vars = c("beta","sigma"))
    row.names(summ)[substr(row.names(summ),1,1) == "b"] <- c("Intercept",x_names)
    row.names(summ)[row.names(summ) == "sigma"] <- "Individual Variance"
  }else{
    summ <- summary(fit, vars = c("beta"))
    row.names(summ) <- c("Intercept",x_names)
  }

  results <- list(summary = summ, distances = d, z = z, jags_model = fit, response = response, effort = effort)

  class(results) <- "latsoc"

  cat("Finished!")

  return(results)

}


#' Plotting latent space models
#'
#' Plots the fitted network object with node positions determined by latent space modelling. 1 and 2 dimensional fits can be plotted directly, while > 3 dimensional fits use multidimensional scaling to project distances into 2 dimensions.
#'
#' @param object A model returned by \code{latent_space}
#' @param post.method Character specifying the method for estimating the posterior latent positions. One of "mean", "median", or "mds", the latter of which uses multidimensional scaling.
#' @param edge.col Square matrix indicating the colors to plot edges. If NULL, edges are plotted in black, with transparency determined by their edge weight.
#' @param edge.lwd Square matrix indicating edge line widths. If NULL, edge widths are proportional to edge weight.
#' @param vertex.pch Plotting characters for vertices/nodes.
#' @param vertex.cex Expansion factor for vertices/nodes.
#' @param vertex.col A vector indicating the colors for plotting nodes.
#' @param labels Logical, should node labels be plotted?
#' @param label.col Colors for vertex labels.
#' @param label.cex Numeric, character expansion factor for labels.
#' @param ... Additional arguments to be passed to \code{plot}
#'
#' @details This function uses base R plotting methods. More advanced network plotting methods are available (e.g. through the ggraph package), but this interface provides an easy way to visualize the results of the model quickly.
#'
#' @export

plot.latsoc <- function(object, post.method = "mean",
                        edge.col = NULL, edge.lwd = NULL,
                        vertex.pch = 21, vertex.col = "black", vertex.cex = 1,
                        labels = F, label.col = NULL, label.cex = NULL,
                        xlim = NULL, ylim = NULL){

  n <- nrow(object$response)

  if(!post.method %in% c("mean", "median", "mds")){
    cat("Invalid post.method, using mean")
    cat("\n")
    post.method <- "mean"
  }

  if(dim(object$z)[[3]] > 2){
    post.method <- "mds"
  }

  if(post.method == "mean"){
    z <- apply(object$z,c(2,3),mean)
  }

  if(post.method == "median"){
    z <- apply(object$z,c(2,3),median)
  }

  if(post.method == "mds"){
    d <- apply(object$distances, c(2,3), mean)
    z <- MASS::isoMDS(as.dist(d))$points
  }

  rate <- object$response/object$effort
  diag(rate) <- 0

  if(is.null(edge.col)){
    edge.col <- rate - min(rate[rate > 0],na.rm=T) + 0.001
    edge.col <- edge.col/max(edge.col,na.rm=T)
    edge.col[edge.col < 0] <- 0
    edge.col <- matrix(rgb(0,0,0,edge.col),nrow=n,ncol=n)
  }

  if(is.null(edge.lwd)){
    edge.lwd <- rate - min(rate[rate > 0],na.rm=T) + 0.001
    edge.lwd <- edge.lwd/max(edge.lwd,na.rm=T)
    edge.lwd <- edge.lwd*4
  }

  plot(z, xlab = "Dimension 1", ylab = "Dimension 2", cex = vertex.cex, bg = vertex.col, pch = 21, ylim = ylim, xlim = xlim)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(rate[i,j] > 0){
        lines(x = z[c(i,j),1], y = z[c(i,j),2], lwd = edge.lwd[i,j], col = edge.col[i,j])
      }
    }
  }
  if(labels){
    names <- colnames(object$response)
    if(!is.null(names)){
      text(x = z[,1], y = z[,2], labels = names, col = label.col, cex = label.cex)
    }
  }

  points(z, pch = 21, cex = vertex.cex, bg = vertex.col)

}
