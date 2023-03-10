% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PoARXEstimation.R
\name{marginalPoARXLogLikelihood}
\alias{marginalPoARXLogLikelihood}
\alias{marginalPoARXScore}
\alias{marginalPoARXInfoMat}
\title{Marginal PoARX log-likelihood computations}
\usage{
marginalPoARXLogLikelihood(theta, obs, exogenous = NULL,
  intercept = TRUE, ylags = NULL, mulags = NULL, yinit = NULL,
  muinit = NULL, ...)

marginalPoARXScore(theta, obs, exogenous = NULL, intercept = TRUE,
  ylags = NULL, mulags = NULL, yinit = NULL, muinit = NULL,
  sameParams = FALSE, ...)

marginalPoARXInfoMat(filter, exogenous, intercept = TRUE, ylags = NULL,
  mulags = NULL, ...)
}
\arguments{
\item{theta}{A vector of parameter values.}

\item{obs}{A matrix of discrete-valued time series observations. Each
column of the matrix represents a separate time series.}

\item{exogenous}{A 3-dimensional array containing the values of the exogenous
covariates. Each new value in dimension 1 represents another observation
of the time series. Each new vector in dimension 2 corresponds to a
new covariate and each matrix of the third dimension represents a
separate time series.}

\item{intercept}{A logical indicator for whether an intercept should be
included in the model.}

\item{ylags}{A vector of time lags in the filter equation relating to the
observations.}

\item{mulags}{A vector of time lags in the filter equation relating to the
intensities.}

\item{yinit}{A matrix detailing the initialisations of the observations.}

\item{muinit}{A matrix detailing the initialisations of the means.}

\item{...}{currently not used.}

\item{sameParams}{A logical argument indicating whether the time series
should be modelled using the same set of PoARX parameter values.}

\item{filter}{A list containing observations and mean values, output from
\code{PoARXFilter}.}
}
\value{
for \code{marginalPoARXLogLikelihood}, a numeric value representing
    the log-likelihood of the observations given the parameters.

for \code{marginalPoARXScore}, a numeric vector representing the
    score function of the observations given the parameters.

for \code{marginalPoARXInfoMat},
   \strong{TODO: Georgi deleted the wrong spec that was here.}
}
\description{
Marginal log-likelihood computations for PoARX models.
}
\details{
\code{marginalPoARXLogLikelihood} computes the log-likelihood for
    observations from the marginal distributions of PoARX models.

\code{marginalPoARXScore} computes marginal PoARX score function,
    the score function for observations from the marginal distributions of
    PoARX models.

\code{marginalPoARXInfoMat} computes the marginal PoARX
    information matrix, the information matrix for observations from the
    marginal distributions of PoARX models.
}
