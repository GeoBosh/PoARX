#' @title Filter mean values for time series following a PoARX model
#'
#' @description
#' Using the observations and initialisation values, the PoARX filter
#' equation estimates the unobserved mean values.
#'
#' @param theta A vector of parameter values.
#' @param obs A matrix of discrete-valued time series observations. Each
#'    column of the matrix represents a separate time series.
#' @param exogenous A 3-dimensional array containing the values of the exogenous
#'    covariates. Each new value in dimension 1 represents another observation
#'    of the time series. Each new vector in dimension 2 corresponds to a
#'    new covariate and each matrix of the third dimension represents a
#'    separate time series.
#' @param intercept A logical indicator for the presence of an intercept in the model.
#' @param ylags A vector of time lags in the filter equation relating to the
#'    observations.
#' @param mulags A vector of time lags in the filter equation relating to the
#'    intensities.
#' @param yinit A matrix detailing the initialisations of the observations.
#' @param muinit A matrix detailing the initialisations of the means.
#' @param ... currently not used.
#'
#' @return A list containing the following:
#'     \item{Y}{A list of three elements: 'All', 'Initial', and 'Observed'.
#'        Each contains a matrix with the relevant y-values.}
#'     \item{Lambda}{A list of three elements: 'All', 'Initial', and 'Calculated'.
#'        Each contains a matrix with the relevant mean values.}
#'
#' @description
#' The PoARX filter equation is
#'    \deqn{\lambda_t = \omega + \sum \alpha_l Y_{t-j_l} + \sum \beta_l
#'       \lambda_{t-j_l} + \eta x_{t-1}}
#'    where \eqn{\lambda} represents the mean of the process, \eqn{Y} the
#'    observations and \eqn{x} the covariates. \eqn{\omega}, {\eqn{\alpha_i}},
#'    {\eqn{\beta_i}} and \eqn{\eta} are all positive parameters that must be
#'    supplied.
#'
#' @export
PoARXFilter <- function(theta, obs, exogenous = NULL,
                        intercept = TRUE, ylags = NULL, mulags = NULL,
                        yinit = NULL, muinit = NULL, ...){
  ## Perform checks
  if(!is.vector(theta))
    stop("'theta' must be a vector")
  if(!is.matrix(obs))
    stop("'obs' must be a matrix")

  if(is.data.frame(yinit))
    yinit <- as.matrix(yinit)
  if(is.data.frame(muinit))
    muinit <- as.matrix(muinit)

  if(!is.null(yinit) & !is.matrix(yinit))
    stop("'yinit' must be a matrix or data.frame")
  if(!is.null(muinit) & !is.matrix(muinit))
    stop("'muinit' must be a matrix or data.frame")

  K <- ncol(obs)
  n <- nrow(obs)

  if(any(dim(exogenous)[1] != n, !dim(exogenous)[3] %in% c(0,K)))
    stop("dimension of 'obs' and 'exogenous' do not match")

  w <- as.numeric(intercept)
  P <- length(ylags)
  Q <- length(mulags)
  r <- max(dim(exogenous)[2], 0)
  if(length(theta) != (w + P + Q + r))
    stop("theta not of correct dimension")

  pmax <- max(ylags, 0)
  qmax <- max(mulags, 0)
  if(is.null(yinit) & pmax > 0)
    stop("'yinit' must be provided")
  if(is.null(muinit) & qmax > 0)
    stop("'muinit' must be provided")

  if(all(P > 0, nrow(yinit) != pmax))
    stop("'yinit' has incorrect number of rows")
  if(all(Q > 0, nrow(muinit) != qmax))
    stop("'muinit' has incorrect number of rows")

  ## Filter the means
  P1 <- min(P,1)
  Q1 <- min(Q,1)
  r1 <- min(r,1)
  Y <- rbind(yinit, obs)

  ## create 'lam' before the loop and don't use rbind() to each row in the loop
  lam <- matrix(0, nrow = 1, ncol = K)
  if(!is.null(muinit)){
    lam <- rbind(lam, muinit)
  }

  if(P > 0)
    thetaP <- theta[(w+P1):(w+P)]
  if(Q > 0)
    thetaQ <- theta[(w+P+Q1):(w+P+Q)]
  if(r > 0)
    thetar <- theta[(w+P+Q+r1):(w+P+Q+r)]

  currow <- nrow(lam)
  lam <- rbind(lam, matrix(0, nrow = n, ncol = K))
  pmqy <- pmax - qmax - ylags  # use below t + pmqy = t + pmax - qmax - ylags
  for(t in (seq_len(n) + qmax)){
    ## t+pmax-qmax-ylags = (t-qmax) + (pmax-ylags) > 0, so no need for pmax
    ## t-mulags > t - rep(qmax, length(mulags) > 0    ,     "
    ## t-qmax > 0
    currow <- currow + 1
    ## changing to if, since otherwise the second terms in '*' are evaluated, eg:
    ##     0 * {print("Argh"); 3}
    ## (potential for trouble and slowing down calc.)
    lam[currow, ] <-
      # theta[1] * w +
      (if(w > 0) theta[1] else 0) +
      (if(P > 0) .colSums(thetaP * Y[t + pmqy, , drop = FALSE], P, K) else 0) +
      (if(Q > 0) .colSums(thetaQ * lam[t - mulags, , drop = FALSE], Q, K) else 0) +
      # TODO: check that r (number of rows) is correct here:
      (if(r > 0) .colSums(matrix(thetar * exogenous[t - qmax, , ],
                                 ncol = K, byrow = TRUE), r, K) else 0)
    ## lam <- rbind(lam, newlam)
  }

  rownames(lam) <- NULL
  lam <- lam[-1 , , drop = FALSE]

  ## Return all observations and mean values
  list(Y = list(All = Y,
                Initial = yinit,
                Observed = obs),
       Lambda = list(All = lam,
                     Initial = muinit,
                     Calculated = lam[(qmax + 1):nrow(lam) , , drop = FALSE]))
}
