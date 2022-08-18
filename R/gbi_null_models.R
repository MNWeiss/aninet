#' Generate Null Distributions for Group-by-Individual Data
#'
#' Produce null distributions of target statistics by permuting a group-by-individual matrix using the MCMC routine proposed by Bejder et al. (1998).
#'
#' @param data A group by individual matrix
#' @param ind_constraint A vector of individual characteristics to constrain permutations; defaults to NULL
#' @param group_constraint A vector group characteristics (usually time windows) to constrain permutations; defaults to NULL
#' @param samples The number of samples to draw for each chain. Defaults to 1000.
#' @param thin The number of steps between each draw. Defaults to 100.
#' @param burnin The number of samples before target statistics are calculated. Defaults to 1000.
#' @param chains The number of independent chains to run. Defaults to 2.
#' @param FUN A function taking a GBI as an input and outputting a vector of any length. See Details.
#' @param ... Other arguments to be passed to FUN
#'
#' @details For some applications, it may be desirable to compare the observed social structure to a null model in which there are no social preferences. This is particularly useful when investigating whether individuals have social preferences, but may also have other uses.
#' In these cases, the canonical method for generating the null model is a Markov Chain Monte Carlo procedure. For each iteration, a change in the matrix is proposed; if the change maintains row an column totals, and falls within the specified constraints, it is accepted, otherwise the matrix stays at its current state.
#' While this method is widely used in social network analysis, and has been implemented in other software, it is often not implemented properly to generate a valid MCMC chain.
#' In addition, none of the current implementations utilize the toolkit that has been developed to check MCMC chain convergence, or provide confidence intervals for p-values. This implementation fills in these gaps.
#'
#' This function takes an argument FUN, which allows users to define custom test statistics. This argument should be a function taking a single argument (the group by individual matrix) and outputting a vector (the test statistic(s)).
#' By default, this function will calculate the test statistics suggested by Whitehead (2008) for testing for non-random social structure: the mean, SD, CV, and portion non-zero association indices (the SRI by default).
#' Other network statistics can be calculated as well (such as clustering coefficient, modularity, etc.).
#'
#' @return An object of class \code{gbi_null}
#'
#' @export

gbi_MCMC <- function(data,
                     ind_constraint = NULL,
                     group_constraint = NULL,
                     samples = 1000,
                     thin = 100,
                     burnin = 1000,
                     chains = 2,
                     FUN = NULL,
                     ...){

  if(!is.matrix(data) | any(!c(data) %in% c(0,1)) | any(is.na(data))){
    stop("Data must be a matrix containing only 1s and 0s")
  }

  N <- ncol(data)
  G <- nrow(data)

  if(is.null(group_constraint)) group_constraint <- rep(1,G)
  if(is.null(ind_constraint)) ind_constraint <- rep(1,N)

  if(!is.vector(group_constraint) | length(group_constraint) != G | any(is.na(group_constraint))) stop("group_constraint must be a vector with length equal to the number of groups")
  if(!is.vector(ind_constraint) | length(ind_constraint) != N | any(is.na(ind_constraint))) stop("ind_constraint must be a vector with length equal to the number of individuals")

  if(is.null(FUN)){
    FUN <- function(gbi){

      x <- get_numerator(gbi,data_format="GBI",return="vector")
      d <- get_denominator(gbi,data_format="GBI",return="vector")
      sri <- x/d

      res <- c(mean(sri,na.rm=T),stats::sd(sri,na.rm=T),stats::sd(sri,na.rm=T)/mean(sri,na.rm=T),mean(sri > 0, na.rm=T))
      names(res) <- c("Mean","SD","CV","Non-zero")
      return(res)

    }
  }
  if(!is.function(FUN)) stop("FUN must be a function")

  min_iters <- gbi_MCMC_iters(data, target_samples=1000, quiet=TRUE)
  if (samples * thin < min_iters) {
    message("Warning: Insufficient iterations for effective sample size of 1000. Minimum number of iterations required is ",
            formatC(min_iters, format="e", digits=1), ". You are only using ", samples, " x ", thin, " = ",
            formatC(samples * thin, format="e", digits=1), " iterations. Null samples will not be reliable.")
  }

  observed <- FUN(data,...)

  chain_res <- list()

  final_gbis <- list()

  for(k in 1:chains){

    cat(paste("Running Chain", k), sep = " ")
    cat("\n")

    gbi.p <- data

    res_matrix <- matrix(nrow = samples, ncol = length(observed))
    colnames(res_matrix) <- names(observed)

    cat(paste("Burning in for", burnin, "iterations", sep = " "))
    cat("\n")

    pb <- utils::txtProgressBar(min = 0, max = burnin*thin,style=3)
    for(i in 1:(burnin*thin)){
      utils::setTxtProgressBar(pb,i)
      cols <- sample(N,2)
      rows <- sample(G,2)
      trial_matrix <- gbi.p[rows,cols]
      if( all(rowSums(trial_matrix) == 1) &
          all(colSums(trial_matrix) == 1) &
          group_constraint[rows[1]] == group_constraint[rows[2]] &
          ind_constraint[cols[1]] == ind_constraint[cols[2]]){

        trial_matrix <- ifelse(trial_matrix == 1, 0, 1)
        gbi.p[rows,cols] <- trial_matrix

      }
    }
    close(pb)

    cat(paste("Sampling for",samples,"iterations",sep=" "))
    cat("\n")

    pb <- utils::txtProgressBar(min = 0, max = samples, style=3)

    for(i in 1:samples){
      utils::setTxtProgressBar(pb,i)
      for(j in 1:thin){
        cols <- sample(N,2)
        rows <- sample(G,2)
        trial_matrix <- gbi.p[rows,cols]
        if( all(rowSums(trial_matrix) == 1) &
            all(colSums(trial_matrix) == 1) &
            group_constraint[rows[1]] == group_constraint[rows[2]] &
            ind_constraint[cols[1]] == ind_constraint[cols[2]]){

          trial_matrix <- ifelse(trial_matrix == 1, 0, 1)
          gbi.p[rows,cols] <- trial_matrix

        }
      }
      res_matrix[i,] <- FUN(gbi.p,...)
    }
    close(pb)
    cat("\n")

    final_gbis[[k]] <- gbi.p
    chain_res[[k]] <- coda::mcmc(res_matrix, thin = thin, start =(burnin+1)*thin, end = thin*(samples+burnin))

  }

  chain_res <- coda::mcmc.list(chain_res)

  results <- list(
    FUN = FUN,
    ind_constraint = ind_constraint,
    group_constraint = group_constraint,
    control = list(thin = thin, samples = samples, burnin =burnin),
    permuted_data = final_gbis,
    observed = observed,
    mcmc = chain_res
  )

  class(results) <- "gbi_null"

  return(results)

}

#' Extend a MCMC Run for Group-By-Individual Permutations
#'
#' Runs additional permutation steps from a \code{gbi_MCMC} run.
#'
#' @param x A \code{gbi_null} object as generated by \code{gbi_MCMC}
#' @param samples The number of additional samples to draw
#' @param ... Additional arguments to generate the test statistic
#'
#' @details This function will extend the sampling of a \code{gbi_null} object. All settings, including the thin and statistics being calculated, will be the same as the initial run.
#'
#' @export

extend.gbi_MCMC <- function(x,samples,...){

  if(samples %% 1 != 0) stop("Argument samples must be an integer")
  if(!is(x,"gbi_null")) stop("Can only extend gbi_null objects")

  observed <- x$observed

  initial_samples <- x$control$samples

  chains <- length(x$mcmc)

  burnin <- x$control$burnin
  thin <- x$control$thin

  N <- ncol(x$permuted_data[[1]])
  G <- nrow(x$permuted_data[[1]])

  group_constraint <- x$group_constraint
  ind_constraint <- x$ind_constraint

  chain_res <- list()
  final_gbis <- list()

  FUN <- x$FUN

  for(k in 1:chains){

    cat(paste("Running Chain", k), sep = " ")
    cat("\n")

    gbi.p <- x$permuted_data[[k]]

    res_matrix <- matrix(nrow = samples, ncol = length(observed))
    colnames(res_matrix) <- names(observed)

    cat(paste("Sampling for",samples,"iterations",sep=" "))
    cat("\n")

    pb <- utils::txtProgressBar(min = 0, max = samples, style=3)

    for(i in 1:samples){
      utils::setTxtProgressBar(pb,i)
      for(j in 1:thin){
        cols <- sample(N,2)
        rows <- sample(G,2)
        trial_matrix <- gbi.p[rows,cols]
        if( all(rowSums(trial_matrix) == 1) &
            all(colSums(trial_matrix) == 1) &
            group_constraint[rows[1]] == group_constraint[rows[2]] &
            ind_constraint[cols[1]] == ind_constraint[cols[2]]){

          trial_matrix <- ifelse(trial_matrix == 1, 0, 1)
          gbi.p[rows,cols] <- trial_matrix

        }
      }
      res_matrix[i,] <- FUN(gbi.p, ...)
    }
    close(pb)
    cat("\n")

    res_matrix <- rbind(as.matrix(x$mcmc[[k]]),res_matrix)

    final_gbis[[k]] <- gbi.p
    chain_res[[k]] <- coda::mcmc(res_matrix, thin = thin, start =(burnin+1)*thin, end = thin*(initial_samples+samples+burnin))

  }

  x$mcmc <- coda::mcmc.list(chain_res)
  x$permuted_data <- final_gbis
  x$control$samples <- x$control$samples + samples

  return(x)

}

#' Plot Null Distributions for Group-By-Individual Matrices
#'
#' Diagnostic plots for GBI null hypotheses tested with the Manly-Bejder MCMC procedure
#'
#' @param x A \code{gbi_null} object as produced by \code{gbi_MCMC}
#'
#' @details This function will first plot the trace plots for the null distributions of each test statistic. It will then plot the change in estimated p-values over the chains.
#'
#' @export

plot.gbi_null <- function(x){
  plot(x$mcmc)
  graphics::par(mfrow = c(which.min(c(4,length(x$observed))),1))
  for(i in 1:length(x$observed)){
    for(k in 1:length(x$mcmc)){
      pval_gr <- cumsum(x$mcmc[[k]][,i] > x$observed[i])/(1:length(x$mcmc[[k]][,i]))
      if(k == 1){
        plot(pval_gr, main = names(x$observed)[i], type = "l", col = k, xlab = "Iteration", ylab = "P(Random > Observed)", ylim = c(0,1))
      }else{
        graphics::points(pval_gr, type = "l", col = k)
      }
    }
  }
}

#' Summarise results of Group-By-Individual Randomisations
#'
#' Generate relevant summaries of null hypothesis testing using MCMC randomisations of group-by-individual matrices
#'
#' @param x a \code{gbi_null} object as produced by \code{gbi_MCMC}
#'
#' @details The function first produces a summary of the null distribution, with mean, median, SD, confidence intervals, and diagnostics (effective sample size and potential scale reduction factor).
#' The function then shows estimated upper, lower, and two-tailed p-values of the observed values against the null. The output includes confidence intervals for the p-values, which are estimated from the effective sample size using binomial theory.
#'
#' @export

summary.gbi_null <- function(x){

  all_mcmc <- do.call(rbind,x$mcmc)

  psrf <- coda::gelman.diag(x$mcmc,multivariate = F)[[1]][,1]

  ess_mean <- mcmcse::ess(x$mcmc)

  null_means <- colMeans(all_mcmc)
  null_median <- apply(all_mcmc,2,stats::median)
  null_sd <- apply(all_mcmc,2,stats::sd)
  CI <- t(apply(all_mcmc,2,stats::quantile,probs=c(0.025,0.975)))

  dist_summary <- cbind(null_means,null_median,null_sd,CI,ess_mean,psrf)
  colnames(dist_summary) <- c("Mean", "Median", "SD", "95% LCI", "95% UCI","ESS","PSRF")

  cat("Null Distribution", "\n")
  print(dist_summary)
  cat("\n")

  ess_pval <- ess_mean

  p_gr <- sapply(1:length(x$observed), function(z){
    mean(all_mcmc[,z] > x$observed[z])
  })
  p_ls <- sapply(1:length(x$observed), function(z){
    mean(all_mcmc[,z] < x$observed[z])
  })
  p_twotail <- sapply(1:length(x$observed), function(z){
    min(c(p_gr[z],p_ls[z]))*2
  })

  ci_gr <- binom::binom.confint(p_gr*ess_pval, ess_pval, method = "exact")[,5:6]
  ci_ls <- binom::binom.confint(p_ls*ess_pval, ess_pval, method = "exact")[,5:6]
  ci_tt <- binom::binom.confint(p_twotail*ess_pval, ess_pval, method = "exact")[,5:6]

  upper_summary <- cbind(p_gr,ci_gr,ess_pval)
  lower_summary <- cbind(p_ls,ci_ls,ess_pval)
  twotail_summary <- cbind(p_twotail,ci_tt,ess_pval)

  colnames(upper_summary) <- colnames(lower_summary) <- colnames(twotail_summary) <- c("P-Value","95% UCI","95% LCI","ESS")
  row.names(upper_summary) <- row.names(lower_summary) <- row.names(twotail_summary) <- names(x$observed)

  cat("P(Random > Observed)","\n")
  print(upper_summary)
  cat("\n")

  cat("P(Random < Observed)","\n")
  print(lower_summary)
  cat("\n")

  cat("Two-Tailed P-values","\n")
  print(twotail_summary)
  cat("\n")

}

#' Minimum number of MCMC iterations needed for Group-by-Individuals randomisations.
#'
#' Uses an analytically-derived equation to estimate the minimum number of MCMC iterations required to generate a target number of effective samples from the null.
#' This function should only be used as a rough guide before running computationally intensive randomisations with \code{gbi_MCMC}. The true number of MCMC iterations may need to be much higher, so MCMC diagnostics should always be used to verify that samples are reliable.
#'
#' @param gbi A group by individual matrix
#' @param target_samples Target number of effective samples
#'
#' @return The minimum number of MCMC iterations required to generate \code{target_samples} effective samples.
#' @export
gbi_MCMC_iters <- function(gbi, target_samples=1000, quiet=FALSE) {
  p = mean(gbi) # Density
  m = length(gbi) # Number of elements

  min_samples <- ceiling((m * target_samples)/(12 * p * (1 - p)))
  if (!quiet) {
    message("Note: The minimum number of iterations reported here may be far lower than the actual number of samples needed. This function is intended to be used as an approximate guide only.")
  }
  return(min_samples)
}
