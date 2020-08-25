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
#' @param vcv.prior Numeric matrix, prior for the fixed effects. Should be a matrix with row and column number equal to the dimensions, representing a parameterization of the Wishart distribution.
#' @param re.prior Numeric vector, prior for individual random effect. Should be the parameters of a Gamma distribution.
#' @param z0 Optional numeric matrix, initial positions of nodes in the latent space.
#' @param ... Further arguments to be passed to runjags.
#'
#' @details Social networks often exhibit transitivity, where a connection between B and C is likely to be stronger if B and C have strong connections to a third individual A. This tendency is often referred to as "triadic closure."
#' One way to account for this in regression settings is to view nodes as being placed within a D dimensional latent space, with edges partially determined by the euclidean distances between them. As euclidean distances are transitive, this modelling framework inherently includes transitivity.
#' This function fits one of these latent space models to a matrix of associations or interactions. For interactions, a Poisson model should be fit, with effort indicating the sampling time per dyad. For associations, a binomial model should be fit, with effort being indicated by the denominator of the association index.
#' In both cases, the response matrix should be a matrix of integers, indicating the number of dyadic interactions or associations.
#' This function fits the model using Gibbs sampling via JAGS and runjags, and allows the user to define priors for fixed effects, the covariance matrix of the latent positions, and individual random effects (if included).
#'
#' @return A list containing the runjags model, matrix of estimated distances, estimated positions, and model summary.
#'
#' @export
latent_space <- function(formula, family = "poisson", dimensions = 2, ind.RE = T, effort, beta.prior = c(0,1e-4), vcv.prior = diag(x = 1, nrow = 2), re.prior = c(0.1,0.1), z0 = NULL, ...){

  if(!is.integer(dimensions)) dimensions <- round(dimensions)
  if(dimensions == 0) stop("Latent space model must have at least 1 latent dimension")

  if(class(vcv.prior) != "matrix" | nrow(vcv.prior) != dimensions){
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

  fit <- runjags::run.jags(model = filestring, data = data, monitor = params, ...)

  res <- do.call(rbind,fit$mcmc)

  cat("Performing transformations on latent positions...")
  cat("\n")

  z <- array(dim = c(nrow(res), n, dimensions))

  if(is.null(z0)) z0 <- MASS::mvrnorm(n=n,mu=rep(0,dimensions),Sigma=diag(1,nrow=dimensions))

  for(i in 1:nrow(res)){
    zi <- matrix(res[i,substr(colnames(res),1,1) == "z"], nrow = n, ncol = dimensions)
    z[i,,] <- MCMCpack::procrustes(zi,z0,translation=T)$X.new
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

  results <- list(summary = summ, distances = d, z = z, jags_model = fit)

  cat("Finished!")

  return(results)

}
