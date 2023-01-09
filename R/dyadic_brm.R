#' Fit Bayesian Dyadic Regression Models
#'
#' This function fits a dyadic regression model (similar to QAP and MRQAP) using Bayesian inference.

#' @param formula A \code{glm} style formula describing the model to be fit. All elements of this formula should be square matrices of the same size.
#' @param family A description of the response distribution and link function. See \link[brms]{brm} and \link[brms]{brmsfamily} for details.
#' @param sampling_effort A square matrix with entries representing dyadic sampling effort. The exact way that this is measured will depend on the data collection protocol (see Details).
#' @param ... Further arguments to be passed to \code{brm}.
#'
#' @details Dyadic regression refers to any analysis of a network where we're interested in how dyadic predictors (e.g. similarlity in some trait) predicts the strength or presence of edges in the network.
#' In this setting, it's reasonable to expect that there might be some dependencies in the response that basic regression doesn't deal with. Imagine if there are individual-level differences in sociality in a social network. All the edges connecting to a very social individual might therefore be higher than expected. This is often referred to as within-row correlation.
#' Traditionally, this has been dealt with through permutation procedures, but this has its drawbacks. While statistical "significance" can be assessed this way, it doesn't help us estimate effect sizes. Furthermore, more complicated models (or even simple ones, like those with interaction effects) don't really work in these cases.
#'
#' A more elegant and useful solution is to actually model the processes we're worried about. We can do this easily in a Bayesian mixed-effect model framework. Here, we model the within-row correlations as multimembership random effects; each dyad gets effects from both of the individuals involved in the dyad.
#' This opens up all kinds of fun possibilities. We can now use interaction effects, as well as smooth terms and random effects in our dyadic regressions.
#'
#' This function uses \code{brms} to do all the heavy lifting; it is just a wrapper to quickly set up the dataset with minimal fuss. In almost all cases, this should be preferred to a permutation method (e.g. \link[aninet]{glmqap}).
#'
#' If \code{sampling_effort = NULL} (the default), no variation in sampling effort is modelled. This is not recommended unless the sampling for each dyad is in fact identical.
#' For association indices (or any network sampled with presence/absence data), the response matrix should be the number of times the pair was observed in association (or interaction), while the number of times they could have been observed in association (the denominator of the index) should be passed to the \code{sampling_effort} argument. For these models, it is recommended to use logit-linear familieslike \code{family = binomial} or \code{family = beta_binomial}.
#' Similarly, for interaction rate data, the number of interactions should be the response matrix, and the observation time should be passed as \code{sampling_effort}. For these types of data, it is recommended to use log-linear count models, like \code{family = poisson} or \code{family = negbinomial}.
#' Any other type of network data (like composite indices, interaction duration, relationship types, etc.) can be used with appropriate families (see \link[stats]{family} and \link[brms]{brmsfamily} for more information). If a matrix is passed to the \code{sampling_effort} argument while using any of these other families, the sampling effort will be treated as proportional weights for each observation, with the assumption that greater sampling effort implies higher certainty in network measurements.
#'
#' The function will incorporate multi-membership random effects in the model by default. Other random effects can be passed as well (e.g. social groups), but will need to be put into square matrices first such that each dyad belongs to one grouping level. See examples.
#' It is also possible to incorporate smooth terms and interaction effects, exactly as they would be passed in a \code{glm} or \code{gam} style formula.
#'
#' Remember, this is a Bayesian model, which means that we're using priors for all of our parameters. It's important to think about what reasonable priors might look like, rather than just using the defaults.
#'
#' @examples
#' # a simple model with a poisson error structure; iterations set quite low for quicker model fitting.
#' model_1 <- dyadic_brm(srkw_contact ~ srkw_kinship, family = poisson, sampling_effort = srkw_sampling, iter = 500)
#' summary(model_1)
#'
#' # a model including the sum and similarity of individual ages
#' age_sum <- sapply(srkw_attributes$age, function(z) z + srkw_attributes$age )
#' age_diff <- sapply(srkw_attributes$age, function(z) abs(z-srkw_attributes$age))
#'
#' model_2 <- dyadic_brm(srkw_contact ~ srkw_kinship + age_sum + age_diff, family = poisson, sampling_effort = srkw_sampling, iter = 500)
#' summary(model_2)
#'
#' # you can add smooth terms as well
#' model_3 <- dyadic_brm(srkw_contact ~ srkw_kinship + age_sum + s(age_diff), family = poisson, sampling_effort = srkw_sampling, iter = 500)
#' summary(model_3)
#'
#' # you can also interact different matrices; for example, if you want to see if there's an interaction between sex and age similarity
#' same_sex <- sapply(srkw_attributes$sex, function(z) ifelse(srkw_attributes$sex == z, 1, 0) )
#'
#' model_4 <- dyadic_brm(srkw_contact ~ srkw_kinship + age_diff*same_sex, family = poisson, sampling_effort = srkw_sampling, iter = 500)
#' summary(model_4)
#'
#' # a model with group-level random effects, using multi-membership terms.
#' mat_1 <- sapply(srkw_attributes$matriline, function(z) rep(z, length(srkw_attributes$matriline)))
#' mat_2 <- t(mat_1)
#' model_5 <- dyadic_brm(srkw_contact ~ srkw_kinship + (1|mm(mat_1,mat_2)), family = poisson, sampling_effort = srkw_sampling, iter = 500)
#' summary(model_5)
#'
#' @return A \code{brmsfit} object.
#'
#' @export

dyadic_brm <- function(formula, sampling_effort = NULL, family = "binomial", ...){

  var_names <- all.vars(formula)

  if(!"character" %in% class(family) & !"family" %in% class(family) & !"function" %in% class(family)){
    stop("'family' must be a character string, family object, or call to a family function")
  }

  if(length(class(family)) == 1){
    if(class(family) == "function"){
      family <- family()
    }
    if(class(family) == "character"){

      if(exists(family)){
        fam_function <- get(family)
      }else{
        fam_function <- getFromNamespace(family, ns = "brms")
      }

      family <- fam_function()

    }
  }

  fam_name <- family$family

  vars <- list()

  for (i in 1:length(var_names)) {
    vars[[i]] <- get(var_names[[i]])
  }

  if(!isSymmetric(vars[[1]])){
    warning("Response matrix is not symmetric; this function is designed for undirected networks.")
  }

  if(!all(unlist(lapply(vars, function(z) nrow(z) == ncol(z) )))){
    stop("Predictors and response must all be square matrices")
  }

  N <- lapply(vars, nrow)
  if(length(unique(N)) != 1) stop("Predictors and response matrices must all have same dimensions")

  if(!is.null(sampling_effort)){
    if(!is.matrix(sampling_effort) | nrow(sampling_effort) != ncol(sampling_effort) | nrow(sampling_effort) != nrow(vars[[1]])){
      stop("Weights must be a square matrix with same dimensions as response")
    }
  }

  id1 <- sapply(1:nrow(vars[[1]]), function(z){
    rep(z,nrow(vars[[1]]))
  })
  id2 <- t(id1)

  df <- do.call(cbind.data.frame,lapply(vars, function(z) z[lower.tri(z)])) #matrix of predictors
  colnames(df) <- var_names
  df$id1 <- id1[lower.tri(id1)]
  df$id2 <- id2[lower.tri(id2)]

  if(!is.null(sampling_effort)){
    df$sampling_effort <- sampling_effort[lower.tri(sampling_effort)]
  }

  formula <- update(formula, .~. + (1|mm(id1,id2)))

  if(!is.null(sampling_effort)){
    if(fam_name %in% c("poisson","negbinomial","zero_inflated_poisson", "zero_inflated_negbinomial","hurdle_poisson", "hurdle_negbinomial")){
      formula <- update(formula, .|rate(sampling_effort) ~ .)
    }else{
      if(fam_name %in% c("binomial", "beta_binomial", "zero_inflated_binomial", "zero_inflated_beta_binomial")){
        formula <- update(formula, .|trials(sampling_effort) ~ .)
      }else{
        formula <- update(formula, .|weights(sampling_effort, scale = TRUE) ~ .)
      }
    }
  }

  bf <- brms::brmsformula(formula)

  model <- brms::brm(bf, data = df, family = family, ...)

  return(model)

}
