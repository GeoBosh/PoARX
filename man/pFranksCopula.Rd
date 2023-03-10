% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CopulaFunctions.R
\name{pFranksCopula}
\alias{pFranksCopula}
\title{Frank's copula CDF for uniform random variables}
\usage{
pFranksCopula(u, rho, ...)
}
\arguments{
\item{u}{a matrix with \eqn{d} columns. Each row specifies a point in the
\eqn{d}-dimensional cube, where the distribution function needs to be
evaluated.}

\item{rho}{a numeric value specifying the dependence parameter.}

\item{...}{currently not used.}
}
\value{
The numeric value(s) of the copula returned as a vector of length
    \code{nrow(u)}.
}
\description{
Frank's copula allows several uniformly distributed random
   variables to be modelled on one probability space.
   % Unclear, so comment out:
   % The copula returns \code{P(U <= u)}.

   This is a wrapper around \code{copula::frankCopula} and related functions
   from package \pkg{copula}.
}
\examples{
#### Generate data
set.seed(20)
u <- matrix(runif(20), ncol = 2)
rho <- 1
pFranksCopula(u, rho)

}
