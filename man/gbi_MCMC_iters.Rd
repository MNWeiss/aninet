% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbi_null_models.R
\name{gbi_MCMC_iters}
\alias{gbi_MCMC_iters}
\title{Minimum number of MCMC iterations needed for Group-by-Individuals randomisations.}
\usage{
gbi_MCMC_iters(gbi, target_samples = 1000, quiet = FALSE)
}
\arguments{
\item{gbi}{A group by individual matrix}

\item{target_samples}{Target number of effective samples}

\item{quiet}{Logical, should messages produced by the function be suppressed?}
}
\value{
The minimum number of MCMC iterations required to generate \code{target_samples} effective samples.
}
\description{
Uses an analytically-derived equation to estimate the minimum number of MCMC iterations required to generate a target number of effective samples from the null.
This function should only be used as a rough guide before running computationally intensive randomisations with \code{gbi_MCMC}. The true number of MCMC iterations may need to be much higher, so MCMC diagnostics should always be used to verify that samples are reliable.
}
