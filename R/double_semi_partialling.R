
dsp.lmer <- function(model, nperm){

  form <- lme4::formula.merMod(model) # get the formula from the model

  if(grepl("*",as.character(form)[[3]], fixed = T) | grepl(":",as.character(form)[[3]], fixed = T)){
    stop("Model contains interaction effects; DSP cannot test these models")
  }

  data <- model@frame # get the data from the model
  data <- cbind(data,rep(1,nrow(data))) # add an intercept column

  colnames(data)[ncol(data)] <- "(Intercept)" # name the intercept column

  fixed <- names(lme4::fixef(model)) #get the names of the fixed effects
  rand <- names(lme4::ranef(model)) #names of the random effects

  X <- data[,fixed,drop=F] # get the fixed predictors as a matrix
  if(any(apply(X,2,class) == "factor")) stop("Fixed predictors must be numeric")
  X <- as.matrix(X) # make the predictors a matrix
  B <- data[,rand,drop=F] #get the random effect matrix

  group.memb <- apply(B,1,function(z){
    paste(z,collapse="-") # collapse all group memberships into a single character vector
  })

  all_groups <- unique(group.memb) # unique group memberships

  if(length(all_groups) == nrow(data)){
    stop("Number of unique random effect levels equal to number of observations")
  }

  Y <- data[,!colnames(data) %in% c(rand,fixed)] #get the responses

  t.perm <- matrix(nrow = nperm, ncol = ncol(X)) #matrix to hold t values

  for(j in 1:ncol(X)){

    x <- X[,j] #get the jth predictor
    z <- X[,-j] #get all other predictors
    resid.x <- stats::residuals(stats::lm(x ~ -1 + z)) #get residuals

    for(i in 1:nperm){ #for each permutation

      for(k in 1:length(all_groups)){ # for each unique group
        resid.x[group.memb == all_groups[k]] <- sample(resid.x[group.memb == all_groups[k]]) #shuffle the residuals in that group
      }

      data.perm <- data #copy original data
      data.perm[,fixed[j]] <- resid.x #plug in residuals

      if(fixed[j] == "(Intercept)"){ # if testing the intercept
        colnames(data.perm)[colnames(data.perm) == fixed[j]] <- "Intercept" # get a nicer intercept
        form.j <- stats::update(form, ~ . + Intercept - 1) # update the formula
        fit.perm <- lme4::lmer(form.j, data = data.perm) # fit the permuted model
        t.perm[i,j] <- summary(fit.perm)$coef["Intercept",3] #save t-value
      }else{
        fit.perm <- lme4::lmer(form, data = data.perm) # fit the permuted model
        t.perm[i,j] <- summary(fit.perm)$coef[j,3] # save the t-value
      }

    }
  }

  t.obs <- summary(model)$coef[,3] #observed t-value

  # calculate p-values
  pval <- sapply(1:length(t.obs), function(z){
    p.gr <- mean(c(t.perm[,z],t.obs[z]) >= t.obs[z])
    p.ls <- mean(c(t.perm[,z],t.obs[z]) <= t.obs[z])
    min(c(p.gr,p.ls))*2
  })

  # make a summary table
  summary_table <- summary(model)$coef
  summary_table <- cbind(summary_table, pval)
  colnames(summary_table)[4] <- "P-value"

  # return the summary table
  return(summary_table)

}

dsp.glmer <- function(model, nperm){

  fam <- lme4::family.merMod(model)
  form <- lme4::formula.merMod(model) # get the formula from the model

  if(grepl("*",as.character(form)[[3]], fixed = T) | grepl(":",as.character(form)[[3]], fixed = T)){
    stop("Model contains interaction effects; DSP cannot test these models")
  }

  data <- model@frame # get the data from the model
  data <- cbind(data,rep(1,nrow(data))) # add an intercept column

  colnames(data)[ncol(data)] <- "(Intercept)" # name the intercept column

  fixed <- names(lme4::fixef(model)) #get the names of the fixed effects
  rand <- names(lme4::ranef(model)) #names of the random effects

  X <- data[,fixed,drop=F] # get the fixed predictors as a matrix
  if(any(apply(X,2,class) == "factor")) stop("Fixed predictors must be numeric")
  X <- as.matrix(X) # make the predictors a matrix
  B <- data[,rand,drop=F] #get the random effect matrix

  group.memb <- apply(B,1,function(z){
    paste(z,collapse="-") # collapse all group memberships into a single character vector
  })

  all_groups <- unique(group.memb) # unique group memberships

  if(length(all_groups) == nrow(data)){
    stop("Number of unique random effect levels equal to number of observations")
  }

  Y <- data[,!colnames(data) %in% c(rand,fixed)] #get the responses

  t.perm <- matrix(nrow = nperm, ncol = ncol(X)) #matrix to hold t values

  for(j in 1:ncol(X)){

    x <- X[,j] #get the jth predictor
    z <- X[,-j] #get all other predictors
    resid.x <- stats::residuals(stats::lm(x ~ -1 + z)) #get residuals

    for(i in 1:nperm){ #for each permutation

      for(k in 1:length(all_groups)){ # for each unique group
        resid.x[group.memb == all_groups[k]] <- sample(resid.x[group.memb == all_groups[k]]) #shuffle the residuals in that group
      }

      data.perm <- data #copy original data
      data.perm[,fixed[j]] <- resid.x #plug in residuals

      if(fixed[j] == "(Intercept)"){ # if testing the intercept
        colnames(data.perm)[colnames(data.perm) == fixed[j]] <- "Intercept" # get a nicer intercept
        form.j <- stats::update(form, ~ . + Intercept - 1) # update the formula
        fit.perm <- lme4::glmer(form.j, data = data.perm, family = fam) # fit the permuted model
        t.perm[i,j] <- summary(fit.perm)$coef["Intercept",3] #save t-value
      }else{
        fit.perm <- lme4::glmer(form, data = data.perm, family = fam) # fit the permuted model
        t.perm[i,j] <- summary(fit.perm)$coef[j,3] # save the t-value
      }

    }
  }

  t.obs <- summary(model)$coef[,3] #observed t-value

  # calculate p-values
  pval <- sapply(1:length(t.obs), function(z){
    p.gr <- mean(c(t.perm[,z],t.obs[z]) >= t.obs[z])
    p.ls <- mean(c(t.perm[,z],t.obs[z]) <= t.obs[z])
    min(c(p.gr,p.ls))*2
  })

  # make a summary table
  summary_table <- summary(model)$coef
  summary_table[,4] <- pval
  colnames(summary_table)[4] <- "P-value"

  # return the summary table
  return(summary_table)

}

dsp.lm <- function(model, nperm){

  form <- stats::formula(model) # get the formula from the model

  if(grepl("*",as.character(form)[[3]], fixed = T) | grepl(":",as.character(form)[[3]], fixed = T)){
    stop("Model contains interaction effects; DSP cannot test these models")
  }

  data <- stats::model.frame(model)
  data <- cbind(data,rep(1,nrow(data))) # add an intercept column

  colnames(data)[ncol(data)] <- "(Intercept)" # name the intercept column

  fixed <- names(stats::coef(model)) #get the names of the fixed effects

  X <- data[,fixed,drop=F] # get the fixed predictors as a matrix
  if(any(apply(X,2,class) == "factor")) stop("Fixed predictors must be numeric")
  X <- as.matrix(X) # make the predictors a matrix

  Y <- data[,!colnames(data) %in% fixed] #get the responses

  t.perm <- matrix(nrow = nperm, ncol = ncol(X)) #matrix to hold t values

  for(j in 1:ncol(X)){

    x <- X[,j] #get the jth predictor
    z <- X[,-j] #get all other predictors
    resid.x <- stats::residuals(stats::lm(x ~ -1 + z)) #get residuals

    for(i in 1:nperm){ #for each permutation

      resid.x <- sample(resid.x) #shuffle the residuals in that group

      data.perm <- data #copy original data
      data.perm[,fixed[j]] <- resid.x #plug in residuals

      if(fixed[j] == "(Intercept)"){ # if testing the intercept
        colnames(data.perm)[colnames(data.perm) == fixed[j]] <- "Intercept" # get a nicer intercept
        form.j <- stats::update(form, ~ . + Intercept - 1) # update the formula
        fit.perm <- stats::lm(form.j, data = data.perm) # fit the permuted model
        t.perm[i,j] <- summary(fit.perm)$coef["Intercept",3] #save t-value
      }else{
        fit.perm <- stats::lm(form, data = data.perm) # fit the permuted model
        t.perm[i,j] <- summary(fit.perm)$coef[j,3] # save the t-value
      }

    }
  }

  t.obs <- summary(model)$coef[,3] #observed t-value

  # calculate p-values
  pval <- sapply(1:length(t.obs), function(z){
    p.gr <- mean(c(t.perm[,z],t.obs[z]) >= t.obs[z])
    p.ls <- mean(c(t.perm[,z],t.obs[z]) <= t.obs[z])
    min(c(p.gr,p.ls))*2
  })

  # make a summary table
  summary_table <- summary(model)$coef
  summary_table[,4] <- pval
  colnames(summary_table)[4] <- "P-value"

  # return the summary table
  return(summary_table)

}

dsp.glm <- function(model, nperm){

  form <- stats::formula(model) # get the formula from the model

  if(grepl("*",as.character(form)[[3]], fixed = T) | grepl(":",as.character(form)[[3]], fixed = T)){
    stop("Model contains interaction effects; DSP cannot test these models")
  }

  data <- stats::model.frame(model)
  data <- cbind(data,rep(1,nrow(data))) # add an intercept column

  colnames(data)[ncol(data)] <- "(Intercept)" # name the intercept column

  fixed <- names(stats::coef(model)) #get the names of the fixed effects

  X <- data[,fixed,drop=F] # get the fixed predictors as a matrix
  if(any(apply(X,2,class) == "factor")) stop("Fixed predictors must be numeric")
  X <- as.matrix(X) # make the predictors a matrix

  Y <- data[,!colnames(data) %in% fixed] #get the responses

  t.perm <- matrix(nrow = nperm, ncol = ncol(X)) #matrix to hold t values

  for(j in 1:ncol(X)){

    x <- X[,j] #get the jth predictor
    z <- X[,-j] #get all other predictors
    resid.x <- stats::residuals(stats::lm(x ~ -1 + z)) #get residuals

    for(i in 1:nperm){ #for each permutation

      resid.x <- sample(resid.x) #shuffle the residuals in that group

      data.perm <- data #copy original data
      data.perm[,fixed[j]] <- resid.x #plug in residuals

      if(fixed[j] == "(Intercept)"){ # if testing the intercept
        colnames(data.perm)[colnames(data.perm) == fixed[j]] <- "Intercept" # get a nicer intercept
        form.j <- stats::update(form, ~ . + Intercept - 1) # update the formula
        fit.perm <- stats::lm(form.j, data = data.perm) # fit the permuted model
        t.perm[i,j] <- summary(fit.perm)$coef["Intercept",3] #save t-value
      }else{
        fit.perm <- stats::lm(form, data = data.perm) # fit the permuted model
        t.perm[i,j] <- summary(fit.perm)$coef[j,3] # save the t-value
      }

    }
  }

  t.obs <- summary(model)$coef[,3] #observed t-value

  # calculate p-values
  pval <- sapply(1:length(t.obs), function(z){
    p.gr <- mean(c(t.perm[,z],t.obs[z]) >= t.obs[z])
    p.ls <- mean(c(t.perm[,z],t.obs[z]) <= t.obs[z])
    min(c(p.gr,p.ls))*2
  })

  # make a summary table
  summary_table <- summary(model)$coef
  summary_table[,4] <- pval
  colnames(summary_table)[4] <- "P-value"

  # return the summary table
  return(summary_table)

}

#' Double-semi-partialling for regression models
#'
#' Perform double-semi-partialling permutations to test the fixed effects of (generalized) linear (mixed effect) models
#'
#' @param model A fitted model object as returned by \code{lm}, \code{glm}, \code{lmer}, or \code{glmer}
#' @param nperm Integer, the number of permutations to perform
#'
#' @details Double-semi-partialling is most commonly used as a permutation procedure for matrix correlations. The basic principles underlying this test, however, are generally applicable to regression models. This is particularly useful for permutation tests with multiple continuous predictors, where simply permuting responses within levels of one or more predictors is not possible.
#' For a given fixed predictor X, we regress X on the remaining fixed predictors Z to arrive at residuals for X. These residuals are then permuted to generate a null distribution. We repeat this for each fixed predictor.
#' For models fit with random effects, this function permutes these predictor residuals within each unique combination of random effects. If the number of unique random effect combinations present is equal to the number of observations, the function returns an error, as no permutation is possible.
#' A drawback of DSP is it does not provide a way to test interaction effects in a principled manner; trying to pass models with interactions to the function will return an error.
#' Note that all categorical variables should be dummy coded into a set of binary variables prior to model fitting.
#' The function currently supports \code{lm} and \code{glm} objects, as well as models fit using \code{lmer} and \code{glmer} in the \code{lme4} package.
#'
#' @return A table with the model estimated coefficients, standard errors, and pivotal statistics, along with the permutation-based P-value.
#'
#' @export

double_semi_partialling <- function(model, nperm){
  cl <- class(model)
  if(length(cl) > 1) cl <- cl[1]
  if(!cl %in% c("glm","lmerMod","glmerMod","lm")){
    stop("Model must be a model fit using lm, glm, lmer, or glmer")
  }
  if(cl == "lmerMod"){
    res <- dsp.lmer(model, nperm)
  }
  if(cl == "glmerMod"){
    res <- dsp.glmer(model, nperm)
  }
  if(cl == "lm"){
    res <- dsp.lm(model,nperm)
  }
  if(cl == "glm"){
    res <- dsp.glm(model,nperm)
  }
  return(res)
}

