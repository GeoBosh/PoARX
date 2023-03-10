% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PoARXPredict.R
\name{PoARXPredict}
\alias{PoARXPredict}
\title{Filter mean values for time series following a PoARX model}
\usage{
PoARXPredict(theta, horizon, nts, obs = NULL, exogenous = NULL,
  intercept = TRUE, ylags = NULL, mulags = NULL, yinit = NULL,
  muinit = NULL, ...)
}
\arguments{
\item{theta}{A vector of parameter values.}

\item{horizon}{A numeric value for the prediction horizon required.}

\item{nts}{A numeric value detailing the number of time series present.}

\item{obs}{An optional matrix of discrete-valued time series observations.
Dimensions are as described in \code{\link{PoARXFilter}}. If these are
not provided and more than one-step ahead is required, they will be
estimated using the intensity values.}

\item{exogenous}{An optional 3-dimensional array containing the values of the
exogenous covariates. Dimensions are as described in
\code{\link{PoARXFilter}}. Must be provided if there are exogenous
covariates in the model.}

\item{intercept}{A logical indicator for the presence of an intercept in the
model.}

\item{ylags}{A vector of time lags in the filter equation relating to the
observations.}

\item{mulags}{A vector of time lags in the filter equation relating to the
intensities.}

\item{yinit}{A matrix detailing the initialisations of the observations.}

\item{muinit}{A matrix detailing the initialisations of the means.}

\item{...}{currently not used.}
}
\value{
TODO: Is this exactly as described in \code{\link{PoARXFilter}}?
}
\description{
Using the observations and the last observed values, use the PoARX filter
    equation to predict the unobserved mean values. Where observations
    are not found, they are estimated using the relevant mean value.

\strong{TODO:} This needs a look. Most of the documentation for this function
    seemed to have been copied verbatim from \code{\link{PoARXFilter}},
    including the title. If this is intended, we can describe them together.
    The descriptions of the parameters below are identical to those in
    \code{\link{PoARXFilter}}, except for \code{obs}.

Description is currently identical to \code{\link{PoARXFilter}}
}
