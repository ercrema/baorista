% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/icarfit.R
\name{icarfit}
\alias{icarfit}
\title{Fits a random walk ICAR model to Aoristic data}
\usage{
icarfit(
  x,
  niter = 1e+05,
  nburnin = 50000,
  thin = 10,
  nchains = 4,
  sigmaPrior = "dexp(rate=1)",
  sigmaSampler = NULL,
  parallel = FALSE,
  seeds = 1:4
)
}
\arguments{
\item{x}{A ProbMat class object}

\item{niter}{Number of MCMC iterations. Default is 500,000.}

\item{nburnin}{Number of iterations discarded for burn-in. Default is 250,000.}

\item{thin}{Thinning interval}

\item{nchains}{Number of MCMC chains}

\item{sigmaPrior}{A string defining prior for the sigma parameter. Default is 'dexp(rate=1)'.}

\item{sigmaSampler}{A list containing settings for the MCMC sampler. Default is null and employs nimble's Default sampler (RW sampler).}

\item{parallel}{Logical specifying whether the chains should be run in parallel or not.}

\item{seeds}{Random seed for each chain. Default is 1:4.}
}
\value{
A \code{fittedICAR} class object containing the original ProbMat class object, posteriors of the probabilities for each time-block and the variance parameter along with their MCMC diagnostics (Gelman Rubin statistic and effective sample size).
}
\description{
Estimates parameters of a multinomial probability distribution that describes the shape of the of the time-frequency distribution of an observed sets of events with chronological uncertainty. The function is wrapper for fitting a 1D random walk ICAR model via nimble.
}
\details{
The function estimates a vector temporally autocorrelated  probability values on the observed data using MCMC as implemented by the nimble package. The model is effectively non-parametric, and at its core it is an implementation of a 1D random ICAR model. Users can specify the prior for the variance parameter through the argument \code{sigmaPrior}. Different settings for this parameter can greatly influence the estimates of the probability vectors. For example using \code{sigmaPrior=dexp(100)} instead of the default \code{sigmaPrior=dexp(1)} would lead to 'flatter' time-series with higher temporal autocorrelation. The distribution parameters defined in \code{sigmaPrior} is also used to generate initialisation values for the MCMC. Please consult the nimble package manual for the syntax required in specifying the prior. MCMC settings such as the choice the sampler, number of iterations, chains, etc can also be specified. Please not that the function is computationally intensive and might require a larger number of iterations to reach satisfactory MCMC convergence.
}
