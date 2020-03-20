#' Estimation of social differentiation
#'
#' Estimates the coefficient of variation of the underlying association probabilities using maximum likelihood, assuming underlying association probabilities follow a beta distribution.
#'
#' @param Num Numeric, numerator of the association indices
#' @param Den Numeric, denominator of association indices
#' @param method Character, indicating what likelihood function to use
#' @param res Resolution used for integration
#' @param initial.params Initial parameters for model fitting
#' @export
#' @details Social differentiation is commonly defined as the coefficient of variation of the true, underlying association probabilities. This estimation procedure assumes that the underlying probabilities follow a beta distribution, and estimates the parameters of this distribution given the observed association indices.
#' The estimation of social differeniation serves as both a useful descriptor of social structure, and a useful measure for determining the power and precision of social analyses. The estimated correlation between true and observed association indices can be derived by dividing the estimated social differentiation by the observed CV of association indices.
#' This function can perform estimation in two ways. The first is using the likelihood function given by Whitehead (2009), using \code{method = "Whitehead"}. The second uses the probability mass function of the beta-binomial distribution using \code{method = "Beta-binomial"}.
#' Because it requires numeric integration, the \code{"Whitehead"} method is much slower than the \code{"Beta-binomial"} method, however has the upside of being directly comparable to results from SOCPROG.
#' Data should be numeric vectors of numerators and denominators, rather than square matrices. See \code{get_numerator()} and \code{get_denominator()} for easy ways to extract these from raw data.
#'
#' @return A matrix containing the estimated social differentiation, the CV of the observed associations, and the estimated correlation between true and observed association indices.
#' @export
social_differentiation <- function(Num, Den, method = c("Whitehead","Beta-binomial"), res = 0.001,  initial.params = c(0.1,0.1)){
  if(length(method) > 1) method <- method[1]
  if(!method %in% c("Whitehead","Beta-binomial")) method = "Whitehead"
  X = Num
  D = Den
  #likelihood functions for social differentiation
  LL.whitehead <- function(z, X, D, delt){
    a = exp(z[1])
    b = exp(z[2])
    deltint <- as.vector(seq(from = delt, to = 1 - delt, by = delt))
    D = as.matrix(D)
    X = as.matrix(X)
    deltint = as.vector(seq(from = delt, to = 1 - delt, by = delt)) #vector of values in [0,1]
    ndel = length(deltint)
    nX = length(X)
    XX = X %*% t(rep(1, ndel))
    DD = D %*% t(rep(1, ndel))
    deldel = rep(1, nX) %*% t(deltint)
    rrx = suppressWarnings(stats::dbeta(deltint, a, b)) #density for beta
    I = ((deldel^XX)*((1-deldel)^(DD-XX))) %*% rrx
    I = sum(log(I))
    -I
  }
  LL.betabinom <- function(z, X, D){
    a = exp(z[1])
    b = exp(z[2])
    ll = VGAM::dbetabinom.ab(X, size = D, shape1 = a, shape2 = b, log = T)
    I = sum(ll)
    -I
  }
  #get parameter estimates
  if(method == "Whitehead"){
    result = stats::optim(initial.params, fn = LL.whitehead, X = X, D = D, delt = res, hessian = T) #MLE for all AIs
  }
  if(method == "Beta-binomial"){
    result = stats::optim(initial.params, fn = LL.betabinom, X = X, D = D, hessian = T) #MLE for all AIs
  }
  a = exp(result$par[1])
  b = exp(result$par[2])
  mean.fit = a/(a+b)
  sd.fit = sqrt((a*b)/((a+b)^2*(a+b+1)))
  estimate = sd.fit/mean.fit
  observed = stats::sd(X/D)/mean(X/D)
  correlation = estimate/observed
  summary = matrix(nrow = 3, ncol = 1)
  row.names(summary) = c("Social Differentiation", "Observed CV", "Correlation")
  summary[1,] = estimate
  summary[2,] = observed
  summary[3,] = correlation
  return(summary)
}

#'Complexity of social associations using mixture models
#'
#'Fits a binomial mixture model to observed associations and calculates the resulting complexity measure.
#'@param Den Numeric, denominator of association indices
#'@param Num Numeric, numerator of association indices
#'@param maxJ Numeric, maximum number of components to fit
#'@param criterion Character, information criteria to use for model selection.
#'@param maxrep Numeric, maximum number of repitions for fitting process
#'@param tolfun Numeric, tolerance for stopping criteria
#'@param minprior Numeric, minimum portion of dyads in a cluster
#'@param maxiter Numeric, maximum number of iterations for convergence
#'
#' @details Social complexity can be defined in a number of ways. This algorithm seeks to measure the diversity of "classes" of assocation occuring at any point. The model describes the association data as a mixture of K association types, each with their own probability of association.
#' By fitting models with K = {1, 2, ..., maxJ}, we can perform model selection to choose the "best" number of relationship types. Each type has a fitted association probability (p) and portion of dyads (a). The portion of association of each type is then q = p*a, and we then measure the entropy of q to arrive at a measure of complexity.
#' The model additionally fits beta-binomial models to derive the dispersion within components.
#' Model selection can use any of AIC, BIC, or ICL to decide on the best model. In simulations ICL performs most consistently, but may underestimate complexity when sample sizes are small.
#'
#' @return A named list with two elements. The \code{summary} element contains a summary of all fit models, including number of components, estimated complexity, and information criteria. The \code{best_model} contains more detailed information on the model chosen by the specified information criterion. This element is a named list, with detailed information on the number of components, component means and frequencies, estimated complexity, likelihood, and information criteria.
#'
association_complexity <- function(Den,Num,maxJ=9,criterion="ICL",maxrep=20, tolfun=1e-6,minprior=0,maxiter=1000){

  bt = function(n,Y,J,maxiter,maxrep,tolfun,minprior){
    fit.rho = function(Q,A,n,Y){
      K = length(Q)
      nY = length(Y)
      nm = matrix(n,ncol=K,nrow=nY)
      Ym = matrix(Y,ncol=K,nrow=nY)
      Am = matrix(A,nrow=nY,ncol=K,byrow=T)
      Qm = matrix(Q,nrow=nY,ncol=K,byrow=T)
      ll.rho = function(r){
        b1 = Qm*(1/r - 1)
        b2 = ((Qm-1)*(r-1))/r
        ll = VGAM::dbetabinom.ab(Ym,nm,b1,b2)*Am
        -sum(log(rowSums(ll)))
      }
      stats::optimise(ll.rho,interval=c(0,1))$minimum
    } #function to fit overdispersion parameter
    mean = list() #lists to hold parameters
    freq = list()
    rho = list()
    mus = list()
    lllq = NULL
    nY = length(Y) #number of dyads
    for(mr in 1:maxrep){
      ll0 = 0
      nm = matrix(n,nrow=nY,ncol=J) #a matrix for the n's
      Ym = matrix(Y,nrow=nY,ncol=J) #a matrix of Y's
      Q = 0.8*stats::runif(J,0,1)
      A = rep(1/J,J)
      K = J
      for(j in 1:maxiter){
        Q = Q[A>minprior] #drop components with no weight (typically not an issue with fuzzy clustering)
        A = A[A>minprior]
        K = length(A) #get new number of components
        nm = matrix(n,ncol=K,nrow=nY)
        Ym = matrix(Y,ncol=K,nrow=nY)
        Am = matrix(A,nrow=nY,ncol=K,byrow=T) #turn component paramters into matrices
        Qm = matrix(Q,nrow=nY,ncol=K,byrow=T)
        ll = stats::dbinom(Ym,nm,Qm)*Am #get likelihoods
        mu = ll/rowSums(ll) #responsibilities
        A = colMeans(mu) #get fractions
        Q = colSums(mu*Y)/colSums(mu*n) #new parameters
        if(any(is.na(Q))|any(is.na(A))) break
        lll = sum(log(rowSums(ll))) #log-likelihood
        if(abs(lll-ll0)<tolfun & j>(maxiter/10)) break #check for convergence
        ll0 = lll
      }
      rho[[mr]] = fit.rho(Q,A,n,Y)
      mean[[mr]] = Q #save parameters
      freq[[mr]] = A
      lllq = c(lllq,lll)
      lllm = max(lllq,na.rm=T)
      if(mr>4 & sum(abs(lllm-lllq)<tolfun,na.rm = T)>4) break
    }
    mean = mean[[which.max(lllq)]]
    freq = freq[[which.max(lllq)]]
    rho = rho[[which.max(lllq)]]
    freq = freq[order(mean)]
    mean = mean[order(mean)]
    K = length(mean)
    Qm = matrix(mean,nrow=nY,ncol=K,byrow=T)
    Am = matrix(freq,nrow=nY,ncol=K,byrow=T)
    Ym = matrix(Y, nrow=nY,ncol=K)
    nm = matrix(n, nrow = nY, ncol= K)
    ll = stats::dbinom(Ym,nm,matrix(mean,nrow=nY,ncol=K,byrow=T))*matrix(freq,nrow=nY,ncol=K,byrow=T)
    mu = ll/rowSums(ll)
    qq = mean*freq
    qq = qq/sum(qq)
    qq = qq[qq>0] #prevents NAs in the entropy estimate
    S = -sum(qq*log(qq))
    lllm = sum(log(rowSums(ll)))
    cl = mu[mu!=0]
    vv = -sum(cl*log(cl))
    AIC = 2*(2*J-1) - 2*lllm
    BIC = log(nY)*(2*J-1) - 2*lllm
    ICL = BIC + 2*vv
    return(list(K.in = J, K.out = length(mean),Mean=mean,Frequency=freq, S = S, rho = rho, logLik=lllm, BIC = BIC, AIC = AIC, ICL = ICL, nrep=mr))
  }

  summary = matrix(ncol = 7,nrow=maxJ) #matrix for summary output
  models = list() #list to hold all models
  colnames(summary) = c("K.in", "K.out", "S", "rho", "AIC","BIC","ICL")
  worse = 0
  for(i in 1:maxJ){
    print(paste("Fitting",i,"Component(s)"))
    res = bt(Den,Num,i,maxrep=maxrep,tolfun=tolfun,minprior=minprior,maxiter=maxiter)
    summary[i,"K.in"] = i
    summary[i,"K.out"] = res$K.out
    summary[i,"S"] = res$S
    summary[i, "rho"] = res$rho
    summary[i,"AIC"] = res$AIC
    summary[i,"BIC"] = res$BIC
    summary[i,"ICL"] = res$ICL
    if(i > 1){
      if(criterion == "AIC"){ #If the current model is worse than the previous, record this.
        worse = ifelse(summary[i,"AIC"]>summary[(i-1),"AIC"], worse+1, 0)
      }
      if(criterion == "BIC"){
        worse = ifelse(summary[i,"BIC"]>summary[(i-1),"BIC"], worse+1, 0)
      }
      if(criterion == "ICL"){
        worse = ifelse(summary[i,"ICL"]>summary[(i-1),"ICL"], worse+1, 0)
      }
    }
    models[[i]] = res
    if(worse > 1) break #If two models in a row have gotten worse, stop fitting
  }
  if(criterion == "BIC"){ #Get the best model, as chosen by the specified criteria
    best = models[[which.min(summary[,"BIC"])]]
  }
  if(criterion == "AIC"){
    best = models[[which.min(summary[,"AIC"])]]
  }
  if(criterion == "ICL"){
    best = models[[which.min(summary[,"ICL"])]]
  }
  summary = summary[1:i,]
  return(list(summary = summary, best_model = best)) #return list containing summary table and best model
}
