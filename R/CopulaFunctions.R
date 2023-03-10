#' @title Frank's copula CDF for uniform random variables
#'
#' @description Frank's copula allows several uniformly distributed random
#'    variables to be modelled on one probability space.
#'    % Unclear, so comment out:
#'    % The copula returns \code{P(U <= u)}.
#'
#'    This is a wrapper around \code{copula::frankCopula} and related functions
#'    from package \pkg{copula}.
#'
#' @param u a matrix with \eqn{d} columns. Each row specifies a point in the
#'     \eqn{d}-dimensional cube, where the distribution function needs to be
#'     evaluated.
#' @param rho a numeric value specifying the dependence parameter.
#' @param ... currently not used.
#'
#' @return The numeric value(s) of the copula returned as a vector of length
#'     \code{nrow(u)}.
#'
#' @examples
#' #### Generate data
#' set.seed(20)
#' u <- matrix(runif(20), ncol = 2)
#' rho <- 1
#' pFranksCopula(u, rho)
#'
#' @export
pFranksCopula <- function(u, rho, ...){
  nr <- nrow(u)
  nc <- ncol(u)
  cop <- frankCopula(param = rho, dim = nc, "TRUE")
  psi(cop, .rowSums(iPsi(cop, u), nr, nc))
}

#' @title Frank's copula PMF for variables from a Poisson model
#'
#' @description Compute the probability mass function for a Frank's copula with
#'     Poisson marginals.
#'
#' @details
#'
#'     Let \eqn{m} be the number of columns of \code{obs}. Each row of
#'     \code{obs} is a vector of observations of a vector r.v. with Poisson
#'     marginals connected by a Frank's copula with parameter \code{rho}.  The
#'     means of the Poisson marginals are in the corresponding row of
#'     \code{means}.
#'
#'     For each row \eqn{o} of \code{obs}, \code{dFranksCopulaPoisson} computes
#'     \eqn{P(Y_1 = o_1, \dots, Y_m = o_m )}, where \eqn{Y_j} are Poisson random
#'     variables with means specified by the corresponding row of \code{means}
#'     and dependence parameter rho.
#'
#'     In this package each column of \code{obs} represents a count time series,
#'     so row \code{t} is a vector of observations at time \code{t}.  The means,
#'     \code{means}, are usually computed by \code{PoARXFilter()}.
#'
#' @param obs a matrix of non-negative integers. Each row represents observed
#'     values of a (vector) count variable, whose probability is to be computed,
#'     see Details.
#' @param means a matrix of non-negative means. It should have the same
#'     dimensions as \code{obs}.
#' @param rho A numeric value specifying the dependence parameter.
#' @param ... currently not used.
#'
#' @return a numeric vector whose i-th element gives the joint probability of
#'     the values in the i-th row of \code{obs}.
#'
#' @seealso \code{\link{pFranksCopula}}, \code{\link{PoARXFilter}}
#' @references \insertRef{Nelsen2006}{PoARX}
#'
#' @export
dFranksCopulaPoisson <- function(obs, means, rho, ...){
  if(any(obs < 0)) stop("Negative observations found")
  if(any(obs %% 1 > 0)) stop("Non-integer observations found")
  ## Ensure 'obs' and 'means' are not data.frames
  if(is.data.frame(obs))
    obs <- as.matrix(obs)
  if(is.data.frame(means))
    means <- as.matrix(means)
  ## Check there's more than one column
  nc <- ncol(obs)
  nr <- nrow(obs)
  if(nc < 2) stop("Cannot form a dependence with less than 2 columns")
  cols <- seq_len(nc)

  ## The joint pmf for n variables consists of 2^n finite differences of the copula
  ## (Nelsen, 1999)
  ## k = l_1 + ... + l_n, eq. 3.17 in Hannah's thesis
  ##
  ## For k = 0 there is only one combination (and it doesn't work in combn)
  ##
  ## 2018-05-17 changed by georgi to create the cop object once.
  ##      (no noticeable effect though)
  cop <- frankCopula(param = rho, dim = nc, "TRUE")

  u <-  ppois(obs, means)
  # FC <- psi(cop, rowSums(iPsi(cop, u))) # was: pFranksCopula(u, rho)
  FC <- psi(cop, .rowSums(iPsi(cop, u), nr, nc))
  for(k in 1:nc){ #### for k > 0, find the relevant combinations of differences
    fc <- combn(cols, k,
                FUN = function(x){
                  obs[ , x] <- obs[ , x] - 1
                  u <- ppois(obs, means)
                  #psi(cop, rowSums(iPsi(cop, u))) # pFranksCopula(u, rho)
                  psi(cop, .rowSums(iPsi(cop, u), nr, nc))
                })
    ## And add/subtract accordingly
    ## @georgi: Number of columns of fc changes depending on k
    FC <- FC + (-1)^k * .rowSums(x = fc, m = nr, n = ncol(fc))
    # FC <- FC + (-1)^k * rowSums(fc)
  }
  ## Remove small negative numbers (machine error)
  FC[FC < 0] <- 0
  FC
}

#' @title Profile log-likelihood w.r.t. the dependence parameter rho
#'
#' @description Profile log-likelihood w.r.t. the dependence parameter rho
#'
#' @details During the IMF process the copula parameters are estimated
#'     separately from the marginal parameters.
#'
#'    \code{rhoFCLogLikelihood} computes the full log-likelihood of Poisson
#'    models coupled using Franks copula with respect to the variation of rho
#'    alone. \strong{(TODO: This sounds bizzare)}
#'
#' @param rho A numeric value specifying the dependence parameter used in
#'    \code{pFranksCopula}.
#' @param obs A matrix of discrete-valued time series observations. Each
#'    column of the matrix represents a separate time series.
#' @param means A matrix of calculated means, from \code{PoARXFilter}. It should
#'    have the same dimensions as \code{obs}.
#' @param zeroFix A numeric value to transform any 0-probability values to.
#' @param ... currently not used.
#'
#' @return for \code{rhoFCLogLikelihood}, a numeric value representing the
#'     log-likelihood of the data for dependence parameter rho.
#'
#' @export
rhoFCLogLikelihood <- function(rho, obs, means, zeroFix = 1e-35, ...){
  probs <- dFranksCopulaPoisson(obs, means, rho)
  probs[probs < zeroFix] <- zeroFix
  sum(log(probs))
}

#' @rdname rhoFCLogLikelihood
#'
#' @details \code{rhoFCScore} computes the profile score function for the
#'     dependence parameter rho, i.e. the derivative of the log-likelihood with
#'     respect to rho.
#'
#' @return for \code{rhoFCScore}, a numeric value representing the score.
#'
#' @export
rhoFCScore <- function(rho, obs, means, zeroFix = 1e-35){
  probs <- dFranksCopulaPoisson(obs, means, rho)
  probsd1 <- .d1drho_dFCPois(obs, means, rho)
  probs[probs < zeroFix] <- zeroFix

  1/probs * probsd1
}

#' @title The profile information matrix for the dependence parameter rho.
#'
#' @description The negative value from twice differentiating the log-likelihood
#'     function for a model coupling Poisson distributions with a Frank's copula
#'     with respect to the dependence parameter.
#'
#' @inheritParams rhoFCLogLikelihood
#'
#' @return A numeric value representing the log-likelihood of the data for
#'    dependence parameter rho.
#'
#' @keywords internal
#' @rdname dot_rhoFCInfoMat
.rhoFCInfoMat <- function(rho, obs, means, zeroFix = 1e-35){
  probs <- dFranksCopulaPoisson(obs, means, rho)
  probsd1 <- .d1drho_dFCPois(obs, means, rho)
  probsd2 <- .d2drho2_dFCPois(obs, means, rho)

  probs[probs < zeroFix] <- zeroFix

  ## left the pmax() for the squared term - probably slightly quicker
  - sum(probsd2/probs - probsd1^2/pmax(probs^2, zeroFix))
}

#' @title Information matrix for the rho-theta part of PoARX model
#'
#' @description
#' Find the information matrix for the theta derivative of the rho derivative
#'    of the PoARX likelihood.
#'
#' @param vartheta A vector of parameter values - all marginal parameters and the
#'    dependence parameter.
#' @param obs A matrix of discrete-valued time series observations. Each
#'    column of the matrix represents a separate time series.
#' @param means A matrix of calculated means, from \code{PoARXFilter}. It should
#'    have the same dimensions as \code{obs}.
#' @param xlist A list containing matrices of exogenous covariates. Each row
#'    represents a new observation, column represents a new covariate. Each element
#'    of the list is for each time series.
#' @param modelSummary A data.frame summarising the features of the PoARX models.
#' @param ylagList A list containing vectors of lags for observations. Each component
#'    of the list represents a time series.
#' @param mulagList A list containing vectors of lags for intensities. Each component
#'    of the list represents a time series.
#' @param yinitList A list of matrices detailing the initialisations of the observations.
#' @param muinitList A list of matrices detailing the initialisations of the means.
#' @param sameParams A logical argument indicating whether all marginal
#'    distributions have the same parameter values.
#' @param zeroFix A numeric value to transform any 0-probability values to.
#' @param ... currently not used.
#'
#' @return A numeric vector representing the score function of the observations
#'    given the parameters.
#'
#' @keywords internal
#' @rdname dot_rhothetaFCInfoMat
.rhothetaFCInfoMat <- function(vartheta, obs, means, xlist = NULL,
                               modelSummary, ylagList = NULL, mulagList = NULL,
                               yinitList = NULL, muinitList = NULL,
                               sameParams = FALSE, zeroFix = 1e-35, ...){
  nr <- nrow(obs)
  nc <- length(vartheta)
  rho <- vartheta[nc]
  L <- dFranksCopulaPoisson(obs, means, rho)
  dLdrho <- .d1drho_dFCPois(obs, means, rho)
  dLdtheta <-.d1dthetaj_dFCPois(vartheta, obs, means, xlist,
                                 modelSummary, ylagList, mulagList,
                                 yinitList, muinitList, sameParams)
  ## dLdtheta <- .dthetaj_dFCPois(vartheta, obs, means, xlist,
  ##                              modelSummary, ylagList, mulagList,
  ##                              yinitList, muinitList, sameParams,
  ##                              type = "orig")
  d2Ldrhotheta <- .d2drhothetaj_dFCPois(vartheta, obs, means, xlist,
                                        modelSummary, ylagList, mulagList,
                                        yinitList, muinitList, sameParams)
  ## d2Ldrhotheta <- .dthetaj_dFCPois(vartheta, obs, means, xlist,
  ##                                  modelSummary, ylagList, mulagList,
  ##                                  yinitList, muinitList, sameParams,
  ##                                  type = "drho")

  ## Change zero probability observations to very small probabilities
  L[L < zeroFix] <- zeroFix
  ## left the pmax() for the squared term - probably slightly quicker
  ## dimension of this is (number of observations) x (number of parameters not rho)
  - .colSums(1/L * d2Ldrhotheta - 1/pmax(L^2, zeroFix) * dLdrho * dLdtheta,
             nr, nc - 1)

}

# first derivative of Frank's copula CDF w.r.t. u_j
.d1duj_pFC <- function(u, rho){
  nr <- nrow(u)
  nc <- ncol(u)
  cop <- frankCopula(param = rho, dim = nc, "TRUE")
  t <- .rowSums(iPsi(cop, u), nr, nc)

  .d1dt_psiFrank(t, rho) * .d1dt_iPsiFrank(u, rho)
}

# first derivative of Frank's copula CDF w.r.t. rho
.d1drho_pFC <- function(u, rho){
  nr <- nrow(u)
  nc <- ncol(u)
  cop <- frankCopula(param = rho, dim = nc, "TRUE")
  t <- .rowSums(iPsi(cop, u), nr, nc)

  res <- .d1dtheta_psiFrank(t, rho) +
    .d1dt_psiFrank(t, rho) *
    .rowSums(.d1dtheta_iPsiFrank(u, rho), nr, nc)

  #### If one of the uniform rvs is 0, then we hit problems (NaN)
  ####   In these cases we should use 0, as the probability in this case is 0.
  res[is.na(res)] <- 0
  res
}

# first derivative of Frank's copula CDF w.r.t. rho and u_j
.d2drhouj_pFC <- function(u, rho){
  nr <- nrow(u)
  nc <- ncol(u)
  cop <- frankCopula(param = rho, dim = nc, "TRUE")
  t <- .rowSums(iPsi(cop, u), nr, nc)

  .d2dthetat_psiFrank(t, rho) * .d1dt_iPsiFrank(u, rho) *
    (1 + .rowSums(.d1dtheta_iPsiFrank(u, rho), nr, nc)) +
    .d1dt_psiFrank(t, rho) * .d2dthetat_iPsiFrank(u, rho)
}

# second derivative of Frank's copula CDF w.r.t. rho (twice)
.d2drho2_pFC <- function(u, rho){
  nc <- ncol(u)
  nr <- nrow(u)
  cop <- frankCopula(param = rho, dim = nc, "TRUE")
  t <- .rowSums(iPsi(cop, u), nr, nc)
  d1rhoiPsi <- .rowSums(.d1dtheta_iPsiFrank(u, rho), nr, nc)

  res <- .d2dtheta2_psiFrank(t, rho) +
    2 * .d2dthetat_psiFrank(t, rho) * d1rhoiPsi +
    .d2dt2_psiFrank(t, rho) * d1rhoiPsi^2 +
    .d1dt_psiFrank(t, rho) *
    .rowSums(.d2dtheta2_iPsiFrank(u, rho), nr, nc)

  #### If one of the uniform rvs is 0, then we hit problems (NaN)
  ####   In these cases we should use 0, as the probability in this case is 0.
  res[is.na(res)] <- 0
  res
}

# first derivative of Frank's copula PMF w.r.t. theta_j for Poisson rvs
.d1dthetaj_dFCPois <- function(vartheta, obs, means, xlist,
                               modelSummary, ylagList, mulagList,
                               yinitList, muinitList, sameParams){
  #### Preliminaries
  # nc <- ncol(obs)     # not needed beyond below
  # cols <- seq_len(nc)
  cols <- seq_len(ncol(obs))
  # N <- length(vartheta) # not needed beyond below
  # rho <- vartheta[N]
  rho <- vartheta[length(vartheta)]

  #### Why do I need newObs and obs?
  ####   We need to keep the original data for past values
  ####   whilst finding rectangle probabilities in the pmf function
  newObs <- obs

  #### See dFranksCopulaPoisson for theory
  FCd1 <- .d1dthetaj_pFCPois(rho, newObs, obs, means, xlist,
                             modelSummary, ylagList, mulagList,
                             yinitList, muinitList, sameParams)
  ## FCd1 <- .dthetaj_pFCPois(rho, newObs, obs, means, xlist,
  ##                            modelSummary, ylagList, mulagList,
  ##                            yinitList, muinitList, sameParams,
  ##                            type = "orig")
  for(k in cols){ #### for k > 0, find the relevant combinations of differences
    fc <- combn(cols, k, FUN = function(x){
      newObs[ , x] <- newObs[ , x] - 1
      .d1dthetaj_pFCPois(rho, newObs, obs, means, xlist,
                         modelSummary, ylagList, mulagList,
                         yinitList, muinitList, sameParams)
      ## .dthetaj_pFCPois(rho, newObs, obs, means, xlist,
      ##                  modelSummary, ylagList, mulagList,
      ##                  yinitList, muinitList, sameParams,
      ##                  type = "orig")
    })
    #### And add/subtract accordingly
    FCd1 <- FCd1 + (-1)^k * apply(fc, 1:2, sum)
  }
  FCd1
}

# first derivative of Frank's copula CDF w.r.t. theta_j for Poisson rvs
.d1dthetaj_pFCPois <- function(rho, newObs, obs, means, xlist,
                               modelSummary, ylagList, mulagList,
                               yinitList, muinitList, sameParams){
  #### Why do I need newObs and obs?
  ####   We need to keep the original data for past values
  ####   whilst finding rectangle probabilities in the pmf function

  #### Preliminary definitions
  nr <- nrow(obs)
  nc <- ncol(obs)
  jseq <- seq_len(nc)
  pmax <- sapply(ylagList, function(x) max(x,0))
  if(length(pmax) == 0)
    pmax <- rep(0, nc)
  qmax <- sapply(mulagList, function(x) max(x,0))
  if(length(qmax) == 0)
    qmax <- rep(0, nc)
  r <- modelSummary$r
  if(nrow(modelSummary) < nc){
    modelSummary[jseq,] <- modelSummary[1,]
  }
  npar <- rowSums(modelSummary)

  #### Create "uniform" variables
  u <-  ppois(newObs, means)

  #### derivative of Frank's copula w.r.t. u_j
  A <- .d1duj_pFC(u, rho)
  A[is.na(A)] <- 0

  #### derivative of u_j w.r.t. lambda_j
  B <- - dpois(newObs, means)

  #### derivative of lambda_j w.r.t. theta_j
  .dlambdajdthetaj <- function(j){
    if(.test_state$get("PoARXFilter") == "0.3.2.9001"){
      C <- NULL
      allObs <- as.matrix(rbind(yinitList[[j]],
                                obs[ , j , drop = FALSE ]))
      allMu <- as.matrix(rbind(muinitList[[j]],
                               means[ , j , drop = FALSE ]))
      if(modelSummary[j,"w"] == 1){
        omega <- 1
      }else{
        omega <- NULL
      }
      for(t in seq_len(nrow(obs))){
        alpha <- allObs[t+pmax[j]-ylagList[[j]],]
        beta <- allMu[t+qmax[j]-mulagList[[j]],]
        eta <- xlist[[j]][t, ]
        derivs <- c(omega, alpha, beta, eta)
        C <- rbind(C, derivs)
      }
      C
    }else{
      allObs <- as.matrix(rbind(yinitList[[j]],
                                obs[ , j , drop = FALSE ]))
      allMu <- as.matrix(rbind(muinitList[[j]],
                               means[ , j , drop = FALSE ]))

      if(modelSummary[j,"w"] == 1){
        omega <- 1
      }else{
        omega <- NULL
      }
      pmy <- pmax[j] - ylagList[[j]]
      qmy <- qmax[j] - mulagList[[j]]

      C <- matrix(0, nrow = nr, ncol = npar[j])
      for(t in seq_len(nr)){
        alpha <- allObs[t + pmy, ]
        beta <- allMu[t + qmy, ]
        eta <- xlist[[j]][t, ]
        C[t, ] <- c(omega, alpha, beta, eta)
      }
      C
    }
  }
  C <- lapply(jseq, .dlambdajdthetaj)

  D <- A * B
  derivs <- lapply(jseq, function(j) D[,j] * C[[j]])
  if(sameParams){
    Reduce("+", derivs)
  }else{
    do.call(cbind, derivs)
  }
}

# second derivative of Frank's copula PMF w.r.t. rho and theta_j for Poisson rvs
.d2drhothetaj_dFCPois <- function(vartheta, obs, means, xlist,
                                  modelSummary, ylagList, mulagList,
                                  yinitList, muinitList, sameParams){
  #### Preliminaries
  # n <- nrow(obs)     # not needed
  # nc <- ncol(obs)    # not needed after next line
  # cols <- seq_len(nc)
  cols <- seq_len(ncol(obs))
  # N <- length(vartheta)    # not needed after next line
  # rho <- vartheta[N]
  rho <- vartheta[length(vartheta)]

  ## Why do I need newObs and obs?
  ##   We need to keep the original data for past values
  ##   whilst finding rectangle probabilities in the pmf function
  newObs <- obs

  ## See dFranksCopulaPoisson for theory
  FCd2 <- .d2drhothetaj_pFCPois(rho, newObs, obs, means, xlist,
                                modelSummary, ylagList, mulagList,
                                yinitList, muinitList, sameParams)
  ## FCd2 <- .dthetaj_pFCPois(rho, newObs, obs, means, xlist,
  ##                          modelSummary, ylagList, mulagList,
  ##                          yinitList, muinitList, sameParams,
  ##                          type = "drho")
  for(k in cols){ #### for k > 0, find the relevant combinations of differences
    fc <- combn(cols, k, FUN = function(x){
      newObs[ , x] <- newObs[ , x] - 1
      .d2drhothetaj_pFCPois(rho, newObs, obs, means, xlist,
                            modelSummary, ylagList, mulagList,
                            yinitList, muinitList, sameParams)
      ## .dthetaj_pFCPois(rho, newObs, obs, means, xlist,
      ##                  modelSummary, ylagList, mulagList,
      ##                  yinitList, muinitList, sameParams,
      ##                  type = "drho")
    })
    #### And add/subtract accordingly
    FCd2 <- FCd2 + (-1)^k * apply(fc, 1:2, sum)
  }
  FCd2
}

# second derivative of Frank's copula CDF w.r.t. rho and theta_j for Poisson rvs
.d2drhothetaj_pFCPois <- function(rho, newObs, obs, means, xlist,
                                  modelSummary, ylagList, mulagList,
                                  yinitList, muinitList, sameParams){
  #### Why do I need newObs and obs?
  ####   We need to keep the original data for past values
  ####   whilst finding rectangle probabilities in the pmf function

  #### Preliminary definitions
  #K <- ncol(newObs)
  nr <- nrow(obs)
  nc <- ncol(obs)
  jseq <- seq_len(nc)

  pmax <- sapply(ylagList, function(x) max(x,0))
  if(length(pmax) == 0)
    pmax <- rep(0, nc)
  qmax <- sapply(mulagList, function(x) max(x,0))
  if(length(qmax) == 0)
    qmax <- rep(0, nc)
  r <- modelSummary$r
  if(nrow(modelSummary) < nc){
    modelSummary[jseq,] <- modelSummary[1,]
  }
  npar <- rowSums(modelSummary)

  #### Create "uniform" variables
  u <-  ppois(newObs, means)

  #### derivative of Frank's copula w.r.t. u_j
  A <- .d2drhouj_pFC(u, rho)
  A[is.na(A)] <- 0

  #### derivative of u_j w.r.t. lambda_j
  B <- - dpois(newObs, means)

  #### derivative of lambda_j w.r.t. theta_j
  ## this is defined in two function environments -
  ##   would it be better to make it global?
  ##   (more arguments would be required)
  .dlambdajdthetaj <- function(j){
    if(.test_state$get("PoARXFilter") == "0.3.2.9001"){
      C <- NULL
      allObs <- as.matrix(rbind(yinitList[[j]],
                                obs[ , j , drop = FALSE ]))
      allMu <- as.matrix(rbind(muinitList[[j]],
                               means[ , j , drop = FALSE ]))
      if(modelSummary[j,"w"] == 1){
        omega <- 1
      }else{
        omega <- NULL
      }
      for(t in seq_len(nrow(obs))){
        alpha <- allObs[t+pmax[j]-ylagList[[j]],]
        beta <- allMu[t+qmax[j]-mulagList[[j]],]
        eta <- xlist[[j]][t, ]
        derivs <- c(omega, alpha, beta, eta)
        C <- rbind(C, derivs)
      }
      C
    }else{
      allObs <- as.matrix(rbind(yinitList[[j]],
                                obs[ , j , drop = FALSE ]))
      allMu <- as.matrix(rbind(muinitList[[j]],
                               means[ , j , drop = FALSE ]))

      if(modelSummary[j,"w"] == 1){
        omega <- 1
      }else{
        omega <- NULL
      }
      pmy <- pmax[j] - ylagList[[j]]
      qmy <- qmax[j] - mulagList[[j]]

      C <- matrix(0, nrow = nr, ncol = npar[j])
      for(t in seq_len(nr)){
        alpha <- allObs[t + pmy, ]
        beta <- allMu[t + qmy, ]
        eta <- xlist[[j]][t, ]
        C[t, ] <- c(omega, alpha, beta, eta)
      }
      C
    }
  }
  C <- lapply(jseq, .dlambdajdthetaj)

  #### Form the derivative of Frank's copula, using chain rule
  D <- A*B
  derivs <- lapply(jseq, function(j) D[,j] * C[[j]])
  if(sameParams){
    Reduce("+", derivs)
  }else{
    do.call(cbind, derivs)
  }
}

## .d1dthetaj_pFCPois and .d2drhothetaj_pFCPois are the same code, but for
##    the definition of the matrix A
##    - create a new function that does both?
##   @georgi: delete as appropriate - this new function or the two similarly
##      coded functions
.dthetaj_pFCPois <- function(rho, newObs, obs, means, xlist,
                             modelSummary, ylagList, mulagList,
                             yinitList, muinitList, sameParams,
                             type = c("orig", "drho")){
  #### Preliminary definitions
  type <- match.arg(type)
  nr <- nrow(obs)
  nc <- ncol(obs)
  jseq <- seq_len(nc)

  pmax <- sapply(ylagList, function(x) max(x,0))
  if(length(pmax) == 0)
    pmax <- rep(0, nc)
  qmax <- sapply(mulagList, function(x) max(x,0))
  if(length(qmax) == 0)
    qmax <- rep(0, nc)
  r <- modelSummary$r
  if(nrow(modelSummary) < nc){
    modelSummary[jseq,] <- modelSummary[1,]
  }
  npar <- rowSums(modelSummary)

  #### Create "uniform" variables
  u <-  ppois(newObs, means)

  #### derivative of Frank's copula w.r.t. u_j
  if(type == "orig"){
    A <- .d1duj_pFC(u, rho)
  }else{
    A <- .d2drhouj_pFC(u, rho)
  }
  A[is.na(A)] <- 0

  #### derivative of u_j w.r.t. lambda_j
  B <- - dpois(newObs, means)

  #### derivative of lambda_j w.r.t. theta_j
  .dlambdajdthetaj <- function(j){
    if(.test_state$get("PoARXFilter") == "0.3.2.9001"){
      C <- NULL
      allObs <- as.matrix(rbind(yinitList[[j]],
                                obs[ , j , drop = FALSE ]))
      allMu <- as.matrix(rbind(muinitList[[j]],
                               means[ , j , drop = FALSE ]))
      if(modelSummary[j,"w"] == 1){
        omega <- 1
      }else{
        omega <- NULL
      }
      for(t in seq_len(nrow(obs))){
        alpha <- allObs[t+pmax[j]-ylagList[[j]],]
        beta <- allMu[t+qmax[j]-mulagList[[j]],]
        eta <- xlist[[j]][t, ]
        derivs <- c(omega, alpha, beta, eta)
        C <- rbind(C, derivs)
      }
      C
    }else{
      allObs <- as.matrix(rbind(yinitList[[j]],
                                obs[ , j , drop = FALSE ]))
      allMu <- as.matrix(rbind(muinitList[[j]],
                               means[ , j , drop = FALSE ]))

      if(modelSummary[j,"w"] == 1){
        omega <- 1
      }else{
        omega <- NULL
      }
      pmy <- pmax[j] - ylagList[[j]]
      qmy <- qmax[j] - mulagList[[j]]

      C <- matrix(0, nrow = nr, ncol = npar[j])
      for(t in seq_len(nr)){
        alpha <- allObs[t + pmy, ]
        beta <- allMu[t + qmy, ]
        eta <- xlist[[j]][t, ]
        C[t, ] <- c(omega, alpha, beta, eta)
      }
      C
    }
  }
  C <- lapply(jseq, .dlambdajdthetaj)

  D <- A * B
  derivs <- lapply(jseq, function(j) D[,j] * C[[j]])
  if(sameParams){
    Reduce("+", derivs)
  }else{
    do.call(cbind, derivs)
  }
}


## .d1dthetaj_dFCPois and .d2drhothetaj_dFCPois are the also the same code, but
##    for FCd1 / FCd2
##    - create a new function that does both?
##   @georgi: delete as appropriate - this new function or the two similarly
##      coded functions
.dthetaj_dFCPois <- function(vartheta, obs, means, xlist,
                             modelSummary, ylagList, mulagList,
                             yinitList, muinitList, sameParams,
                             type = c("orig", "drho")){
  #### Preliminaries
  type <- match.arg(type)
  cols <- seq_len(ncol(obs))
  rho <- vartheta[length(vartheta)]

  #### Why do I need newObs and obs?
  ####   We need to keep the original data for past values
  ####   whilst finding rectangle probabilities in the pmf function
  newObs <- obs

  #### See dFranksCopulaPoisson for theory
  FC <- .dthetaj_pFCPois(rho, newObs, obs, means, xlist,
                         modelSummary, ylagList, mulagList,
                         yinitList, muinitList, sameParams,
                         type)
  for(k in cols){ #### for k > 0, find the relevant combinations of differences
    fc <- combn(cols, k, FUN = function(x){
      newObs[ , x] <- newObs[ , x] - 1
      .dthetaj_pFCPois(rho, newObs, obs, means, xlist,
                       modelSummary, ylagList, mulagList,
                       yinitList, muinitList, sameParams,
                       type)
    })
    #### And add/subtract accordingly
    FC <- FC + (-1)^k * apply(fc, 1:2, sum)
  }
  FC
}

# first derivative of Frank's copula PMF w.r.t. rho for Poisson rvs
.d1drho_dFCPois <- function(obs, means, rho, ...){
  nr <- nrow(obs)
  nc <- ncol(obs)
  if(nc < 2) stop("Cannot form a dependence with less than 2 columns")
  cols <- seq_len(nc)

  #### See dFranksCopulaPoisson for theory
  u <-  ppois(obs, means)
  FCd1 <- .d1drho_pFC(u, rho)
  for(k in cols){ #### for k > 0, find the relevant combinations of differences
    fc <- combn(cols, k, FUN = function(x){
      obs[ , x] <- obs[ , x] - 1
      u <- ppois(obs, means)
      .d1drho_pFC(u, rho)
    })
    #### And add/subtract accordingly
    FCd1 <- FCd1 + (-1)^k * .rowSums(fc, nr, ncol(fc))
  }
  FCd1
}

# second derivative of Frank's copula PMF w.r.t. rho (twice) for Poisson rvs
.d2drho2_dFCPois <- function(obs, means, rho, ...){
  nr <- nrow(obs)
  nc <- ncol(obs)
  if(nc < 2) stop("Cannot form a dependence with less than 2 columns")
  cols <- seq_len(nc)

  #### See dFranksCopulaPoisson for theory
  u <-  ppois(obs, means)
  FCd2 <- .d2drho2_pFC(u, rho)
  for(k in cols){ #### for k > 0, find the relevant combinations of differences
    fc <- combn(cols, k, FUN = function(x){
      obs[ , x] <- obs[ , x] - 1
      u <- ppois(obs, means)
      .d2drho2_pFC(u, rho)
    })
    #### And add/subtract accordingly
    FCd2 <- FCd2 + (-1)^k * .rowSums(fc, nr, ncol(fc))
  }
  FCd2
}
