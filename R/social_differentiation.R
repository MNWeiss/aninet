#' Estimation of social differentiation
#'
#' Estimates the coefficient of variation of the underlying association probabilities using maximum likelihood, assuming underlying association probabilities follow a beta distribution.
#'
#' @param Num Numeric vector, numerator of the association indices
#' @param Den Numeric vector, denominator of association indices
#' @param initial.params Initial parameters for model fitting
#' @param nsim Number of parametric bootstraps to draw for error estimation. See details.
#'
#' @details Social differentiation is commonly defined as the coefficient of variation of the true, underlying association probabilities. This estimation procedure assumes that the underlying probabilities follow a beta distribution, and estimates the parameters of this distribution given the observed association indices.
#' The estimation of social differentiation serves as both a helpful descriptor of social structure, and a useful measure for determining the power and precision of social analyses. The estimated correlation between true and observed association indices can be derived by dividing the estimated social differentiation by the observed CV of association indices.
#' In some edge cases, where social differentiation is very high, a correlation greater than 1 can be estimated. In these cases we recommend primarily interpreting the lower bound of the confidence interval.
#'
#' This function uses the likelihood of the beta-binomial distribution to estimate social differentiation. The original method for social differentiation (and the one used by SOCPROG) uses a slightly different method involving integrating over possible probability values.
#' The beta-binomial method is faster and more precise, but cannot be used when the denominators of the association indices are not integers (e.g. because the HWI was calculated). In these cases, the function will default to the integration method, with a warning.
#'
#' This function will return estimated standard errors and confidence intervals. This is estimated by drawing simulated values using the estimated parameters and their variance-covariance matrix. A more robust estimate of confidence can be achieved throu a bootstrap of the raw data. To save computation time, we recommend setting nsim < 1 in these cases.
#'
#' @return A matrix containing the estimated social differentiation, the CV of the observed associations, and the estimated correlation between true and observed association indices, along with standard errors and confidence intervals.
#'
#' @examples
#' X <- get_numerator(srkw_sightings, return = "vector", data_format = "GBI")
#' D <- get_denominator(srkw_sightings, return = "vector", data_format = "GBI")
#' social_differentiation(X, D)
#'
#' @export
social_differentiation <- function(Num, Den, initial.params = c(0.1,0.1), nsim = 100000){

  X <- Num
  D <- Den

  X <- X[D > 0]
  D <- D[D > 0]

  if(any(X %% 1 != 0)) stop("Numerators contain non-integers")

  #likelihood functions for social differentiation
  LL.betabinom <- function(z, X, D){
    a <- exp(z[1])
    b <- exp(z[2])
    ll <- VGAM::dbetabinom.ab(X, size = D, shape1 = a, shape2 = b, log = T)
    I <- sum(ll)
    -I
  }

  if(any(D %% 1 != 0)){

    warning("Denominators contain non-integers, defaulting to integration method.")

    LL.betabinom <- function(z,X,D){

      a <- exp(z[1])
      b <- exp(z[2])

      int_fun <- function(p,x,d){
        dbeta(p,a,b) * p^x * (1-p)^(d-x) * choose(d, x)
      }

      int <- sapply(1:length(X), function(index){
        integrate(int_fun, lower = 0, upper = 1, x = X[index], d = D[index])$value
      })

      -sum(log(int))

    }

  }

  #get parameter estimates
  result <- stats::optim(initial.params, fn = LL.betabinom, X = X, D = D, hessian = T) #MLE for all AIs
  #transform parameters
  a <- exp(result$par[1])
  b <- exp(result$par[2])
  #calculate mean and standard deviation
  mean.fit <- a/(a+b)
  sd.fit <- sqrt((a*b)/((a+b)^2*(a+b+1)))
  #estimate of social differentiation
  estimate <- sd.fit/mean.fit
  #observed CV
  observed <- stats::sd(X/D)/mean(X/D)
  #estimated correlation
  correlation <- estimate/observed

  if(nsim > 1){
    #sample parameters based on estimates and hessian matrix
    samp <- MASS::mvrnorm(n = nsim, mu = result$par, Sigma = solve(result$hessian))
    #distribution of CVs
    cv.samp <- apply(samp,1,function(z){
      a <- exp(z[1])
      b <- exp(z[2])
      mean.fit = a/(a+b)
      sd.fit = sqrt((a*b)/((a+b)^2*(a+b+1)))
      sd.fit/mean.fit
    })
    #get SEs and CIs
    se_S <- sd(cv.samp, na.rm = T)
    se_r <- sd(cv.samp/observed, na.rm = T)
    ci_S <- quantile(cv.samp, c(0.025,0.975), na.rm = T)
    ci_r <- quantile(cv.samp/observed, c(0.025,0.975), na.rm = T)
  }else{
    se_S <- NA
    se_r <- NA
    ci_S <- NA
    ci_r <- NA
  }

  #make summary table
  summary <- matrix(nrow = 3, ncol = 4)
  row.names(summary) <- c("Observed CV","Social Differentiation", "Correlation")
  colnames(summary) <- c("Estimate", "SE", "Lower CI", "Upper CI")

  summary[1,1] <- observed
  summary[2,1] <- estimate
  summary[2,2] <- se_S
  summary[2,c(3,4)] <- ci_S
  summary[3,1] <- correlation
  summary[3,2] <- se_r
  summary[3,c(3,4)] <- ci_r

  return(summary)
}
