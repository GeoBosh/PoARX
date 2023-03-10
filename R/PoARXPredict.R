#' @title Predict wrapper for the intensity values of a PoARX model
#'
#' @description Used in the \code{predict.PoARX} function and contains the
#'    "hard work" of the predictions. New data can be fed into the model and
#'    a prediction given - this can be done h steps into the future, where
#'    the variance will grow with h. Missing observation values are estimated by
#'    the respective intensity values. If a response is to be returned, there is an
#'    option to form a confidence interval for this value, using simulation.
#'
#' @param object A \code{"PoARX"} object from \code{PoARXFit}.
#' @param newdata An optional data.frame for the new covariates and observations.
#'    The names of the data.frame should match the original variable names. If some new
#'    observations have not been obtained, put an NA value where the missing value is
#'    and the observation will be estimated by \code{PoARXPredict}.
#' @param type A character argument indicating the type of response. If \code{"response"},
#'    return the mean parameter for each time series. If \code{"probability"}, give the
#'    probability of the response.
#' @param interval A logical argument for the calculation of confidence intervals.
#' @param cl A numeric value representing the desired confidence level for predictions of
#'    response values.
#' @param B A number of desired simulations in the calculation of confidence intervals in
#'    any \code{'response'} values.
#' @param trace A logical indicator for the printing of update messages.
#' @param ... currently not used.
#'
#' @return If \code{type = "response"} a data.frame is returned containing
#'    \eqn{\lambda} values for each time series, along with optional confidence intervals.
#'    If \code{type = "probability"} and the respective observations are provided, then
#'    a matrix of probabilities is returned. If no observations are given, there
#'    are too many possible values, so the \eqn{\lambda} values are returned instead,
#'    as if \code{type = "response"}.
#'
#' @keywords internal
#' @rdname dot_predictLambda
.predictLambda <- function(object, newdata){
  #                           interval = FALSE,
  #                           cl = 0.95, B = 1000, trace){
  if(!all(names(object$X) %in% names(newdata))){
    nm1 <- names(object$X)[!names(object$X) %in% names(newdata)]
    stop(paste(nm1, "not found in 'newdata'"))
  }

  #### basic definitions
  w <- object$modelSummary$w
  P <- object$modelSummary$P
  Q <- object$modelSummary$Q
  ylags <- object$ylags   # so we don't keep extracting later
  mulags <- object$mulags # so we don't keep extracting later
  pmax <- max(ylags, 0)
  qmax <- max(mulags, 0)
  r <- object$modelSummary$r
  K <- object$nts
  Kseq <- seq_len(K)
  n1 <- nrow(newdata)
  N <- w + P + Q + r
  #   M <- (object$npar-object$rhoInd)/N

  #### new xlist
  Xnew <- lapply(Kseq, function(j)
    model.part(object$formula, data = newdata, rhs = j))
  Xnew <- lapply(Xnew, function(Xi)
    Xi[,colnames(Xi) != "(Intercept)" , drop = FALSE])

  #### new observation matrix - keep NAs
  Ynew <- NULL
  if(all(names(object$Y) %in% names(newdata))){
    Ynew <- model.part(object$formula, data = newdata, lhs = Kseq)
    Ynew <- as.matrix(Ynew)
  }

  if(.test_state$get("PoARXFilter") == "0.3.2.9001"){
    if(object$sameParams){
      yinit <- object$Y[(object$nobs - pmax + 1):object$nobs,]
      yinit <- as.matrix(yinit)
      muinit <- object$fitted.values[(object$nobs - qmax + 1):object$nobs,]
      muinit <- as.matrix(muinit)

      intercept <- all(as.logical(object$modelSummary$w))
      if(r > 0){
        Xnew <- array(unlist(Xnew), dim = c(n1, r, K))
      }else{
        Xnew <- array(0, dim = c(n1, r, K))
      }

      newlam <- PoARXPredict(object$theta, n1, K, Ynew, Xnew,
                             intercept, ylags[,1], mulags[,1],
                             yinit, muinit)
    }else{
      #### parameters
      mark <- rep(Kseq, times = N)
      theta0 <- object$theta
      if(object$rhoInd)
        theta0 <- theta0[-object$npar]
      thetaj <- split(theta0, mark)
      names(thetaj) <- NULL
      obsinit = NULL
      if(pmax > 0){
        .yinit <- function(j){
          pmax1 <- max(object$ylagList[[j]], 0)
          if(pmax1 > 0){
            yinit <- object$Y[(object$nobs - pmax1 + 1):object$nobs , j , drop = F]
            as.matrix(yinit)
          }else{
            matrix(0, nrow = 0, ncol = 1)
          }

        }
        obsinit <- lapply(Kseq, .yinit)
      }
      meaninit = NULL
      if(qmax > 0){
        .muinit <- function(j){
          qmax1 <- max(object$mulagList[[j]], 0)
          if(qmax1 > 0){
            muinit <- object$fitted.values[(object$nobs - qmax1 + 1):object$nobs , j , drop = F]
            as.matrix(muinit)
          }else{
            matrix(0, nrow = 0, ncol = 1)
          }
        }
        meaninit <- lapply(Kseq, .muinit)
      }
      if(any(nrow(obsinit) == 0))
        obsinit <- NULL
      if(any(nrow(meaninit) == 0))
        meaninit <- NULL

      newlam <- lapply(Kseq, function(j){
        ob <- Ynew[ , j , drop = FALSE]
        r <- object$modelSummary$r[j]
        if(r>0){
          ex <- array(unlist(Xnew[[j]]), dim = c(dim(Xnew[[j]]),1))
        }else{
          ex <- array(0, dim = c(n1,r,1))
        }
        intercept <- as.logical(object$modelSummary$w[j])
        ylags <- object$ylagList[[j]]
        if(length(ylags) == 0)
          ylags <- NULL
        mulags <- object$mulagList[[j]]
        if(length(mulags) == 0)
          mulags <- NULL

        PoARXPredict(theta        =  thetaj[[j]],
                     obs          =  ob,
                     exogenous    =  ex,
                     intercept    =  intercept,
                     ylags        =  ylags,
                     mulags       =  mulags,
                     yinit        =  obsinit[[j]],
                     muinit       =  meaninit[[j]])
      })
      newlam <- do.call(cbind, newlam)
      newlam <- as.data.frame(newlam)
      names(newlam) <- paste("fit", Kseq, sep = ".")
    }
  }else{
    if(object$sameParams){
      yinit <- object$Y[(object$nobs - pmax + 1):object$nobs,]
      yinit <- as.matrix(yinit)
      muinit <- object$fitted.values[(object$nobs - qmax + 1):object$nobs,]
      muinit <- as.matrix(muinit)

      P1 <- min(P,1)
      Q1 <- min(Q,1)
      r1 <- min(r,1)

      Y <- rbind(yinit, Ynew)
      if(nrow(Y) < n1 + pmax)
        Y <- rbind(Y, matrix(NA, nrow = n1 + pmax - nrow(Y), ncol = K))

      lam <- matrix(0, nrow = 1, ncol = K)
      if(!is.null(muinit)){
        lam <- rbind(lam, muinit)
      }

      if(r > 0){
        Xnew <- array(unlist(Xnew), dim = c(n1, r, K))
      }else{
        Xnew <- array(0, dim = c(n1, r, K))
      }

      thetaw <- if(w > 0) object$theta[1] else 0
      if(P > 0)
        thetaP <- object$theta[(w+P1):(w+P)]
      if(Q > 0)
        thetaQ <- object$theta[(w+P+Q1):(w+P+Q)]
      if(r > 0)
        thetar <- object$theta[(w+P+Q+r1):(w+P+Q+r)]

      currow <- nrow(lam)
      lam <- rbind(lam, matrix(0, nrow = n1, ncol = K))
      pmqy <- pmax - qmax - ylags
      for(t in (seq_len(n1) + qmax)){
        currow <- currow + 1
        newlam <- thetaw +
          (if(P > 0) .colSums(thetaP * Y[t + pmqy, , drop = FALSE], P, K)
           else 0) +
          (if(Q > 0) .colSums(thetaQ * lam[t - mulags, , drop = FALSE], Q, K)
           else 0) +
          (if(r > 0) .colSums(matrix(thetar * Xnew[t - qmax, , ],
                                     ncol = K, byrow = TRUE), r, K) else 0)
        lam[currow, ] <- newlam
        naObsInd <- is.na(Y[t-qmax+pmax,])
        if(any(naObsInd))
          Y[t-qmax+pmax,][naObsInd] <- newlam[naObsInd]
      }
      rownames(lam) <- NULL
      # Remove the first row (0s)
      lam <- lam[ -1 , , drop = FALSE]

      #### Return all observations and mean values
      newlam <- lam[(qmax + 1):nrow(lam) , , drop = FALSE]
    }else{
      #### parameters
      mark <- rep(Kseq, times = N)
      theta0 <- object$theta
      if(object$rhoInd)
        theta0 <- theta0[-object$npar]
      thetaj <- split(theta0, mark)
      names(thetaj) <- NULL

      if(pmax > 0){
        .yinit <- function(j){
          pmax1 <- max(object$ylagList[[j]], 0)
          if(pmax1 > 0){
            yinit <- object$Y[(object$nobs - pmax1 + 1):object$nobs , j , drop = F]
            as.matrix(yinit)
          }else{
            NULL
          }
        }
        obsinit <- lapply(Kseq, .yinit)
      }else{
        obsinit <- NULL
      }

      if(qmax > 0){
        .muinit <- function(j){
          qmax1 <- max(object$mulagList[[j]], 0)
          if(qmax1 > 0){
            muinit <- object$fitted.values[(object$nobs - qmax1 + 1):object$nobs , j , drop = F]
            as.matrix(muinit)
          }else{
            NULL
          }
        }
        meaninit <- lapply(Kseq, .muinit)
      }else{
        meaninit <- NULL
      }

      .newlamj <- function(j){
        thetajj <- thetaj[[j]]

        wj <- object$modelSummary$w[j]
        Pj <- object$modelSummary$P[j]
        Qj <- object$modelSummary$Q[j]
        rj <- object$modelSummary$r[j]

        yjlags <- object$ylagList[[j]]
        if(length(yjlags) == 0)
          yjlags <- NULL
        mujlags <- object$mulagList[[j]]
        if(length(mujlags) == 0)
          mujlags <- NULL

        pjmax <- max(yjlags, 0)
        qjmax <- max(mujlags, 0)

        Pj1 <- min(Pj,1)
        Qj1 <- min(Qj,1)
        rj1 <- min(rj,1)

        Yj <- rbind(obsinit[[j]],
                    Ynew[ , j , drop = FALSE])
        if(nrow(Yj) < n1 + pjmax)
          Yj <- rbind(Yj,
                      matrix(NA, nrow = n1 + pjmax - nrow(Yj), ncol = K))

        lamj <- matrix(0, nrow = 1, ncol = 1)
        if(!is.null(meaninit)){
          lamj <- rbind(lamj, meaninit[[j]])
        }

        if(rj > 0){
          ex <- array(unlist(Xnew[[j]]), dim = c(dim(Xnew[[j]]), 1))
        }else{
          ex <- NULL
        }

        thetajw <- if(wj > 0) thetajj[1] else 0
        if(Pj > 0)
          thetajP <- thetajj[(wj+Pj1):(wj+Pj)]
        if(Qj > 0)
          thetajQ <- thetajj[(wj+Pj+Qj1):(wj+Pj+Qj)]
        if(rj > 0)
          thetajr <- thetajj[(wj+Pj+Qj+rj1):(wj+Pj+Qj+rj)]

        currow <- nrow(lamj)
        lamj <- rbind(lamj, matrix(0, nrow = n1, ncol = 1))
        pmqy <- pjmax - qjmax - yjlags
        for(t in (seq_len(n1) + qjmax)){
          currow <- currow + 1
          newlam <- thetajw +
            (if(Pj > 0) .colSums(thetajP * Yj[t + pmqy, , drop = FALSE], Pj, 1)
             else 0) +
            (if(Qj > 0) .colSums(thetajQ * lamj[t - mujlags, , drop = FALSE],
                                Qj, 1) else 0) +
            (if(rj > 0) .colSums(matrix(thetajr * ex[t - qjmax, , ], ncol = 1,
                                       byrow = TRUE), rj, 1) else 0)
          lamj[currow, ] <- newlam
          naObsInd <- is.na(Yj[t - qjmax + pjmax,])
          if(any(naObsInd))
            Yj[t - qmax + pmax,][naObsInd] <- newlam[naObsInd]
        }
        rownames(lamj) <- NULL
        # remove the first row - 0s
        lamj <- lamj[ -1 , , drop = FALSE]

        lamj[(qjmax + 1):nrow(lamj) , , drop = FALSE]
      }
      newlam <- lapply(Kseq, .newlamj)
      newlam <- do.call(cbind, newlam)
      newlam <- as.data.frame(newlam)
      names(newlam) <- paste("fit", Kseq, sep = ".")
    }
  }
  newlam
}

#### I've suspended the interval bit due to comment from Georgi
# .predictLambda <- function(object, newdata, interval = FALSE,
#                            cl = 0.95, B = 1000, trace){
#   #if(!all(names(newdata) %in% names(object$X))){
#   #  nm1 <- names(newdata)[!names(newdata) %in% names(object$X)]
#   #  stop(paste(nm1, "not found in 'PoARX' object"))
#   #}
#   if(!all(names(object$X) %in% names(newdata))){
#     nm1 <- names(object$X)[!names(object$X) %in% names(newdata)]
#     stop(paste(nm1, "not found in 'newdata'"))
#   }
#
#   #### basic definitions
#   P <- length(object$ylags)
#   Q <- length(object$mulags)
#   pmax <- max(object$ylags, 0)
#   qmax <- max(object$mulags, 0)
#   r <- object$nexo
#   K <- object$nts
#   n1 <- nrow(newdata)
#   N <- object$intercept + P + Q + r
#   M <- (object$npar-object$rhoInd)/N
#
#   #### new exogenous array
#   Xdat <- model.part(object$formula, data = newdata, rhs = seq_len(K))
#   Xnew <- array(as.matrix(Xdat), dim = c(n1, r, K))
#
#   #### new observation matrix - keep NAs
#   Ynew <- NULL
#   if(all(names(object$Y) %in% names(newdata))){
#     Ynew <- model.part(object$formula, data = newdata, lhs = seq_len(K))
#     Ynew <- as.matrix(Ynew)
#   }
#
#   #### parameters
#   theta <- object$theta
#   if(object$rhoInd)
#     theta <- theta[-object$npar]
#
#   varCond <- !is.null(object$varCovar)
#   if(!varCond){
#     message("no variance matrix available to simulate")
#   }
#   if(interval & varCond){
#     if(trace)
#       print("intervals require simulation - this may take time")
#
#     coef <- matrix(NA, nrow = B, ncol = (object$npar - object$rhoInd))
#     constr <- PoARXConstraints(M = M, r = r, rho = FALSE,
#                                intercept = object$intercept,
#                                ylags = object$ylags,
#                                mulags = object$mulags)
#     l <- 0
#     while(l < B){
#       coef1 <- mvrnorm(2*(B-l), theta,
#                        object$varCovar[-object$npar, -object$npar])
#       coef1 <- coef1[apply(coef1, 1, function(x)
#         all(constr$ui %*% x >= constr$ci)),]
#       n1 <- nrow(coef1)
#       if(max(n1, 0) > 0){
#         if((B-l) >= n1){
#           coef[(l+1):(l+n1),] <- coef1
#         }else{
#           coef[(l+1):B,] <- coef1[sample(seq_len(nrow(coef1)), B-l),]
#         }
#       }
#       l <- l + max(nrow(coef1), 0)
#     }
#     alpha <- (1-cl)/2
#   }
#
#   #### Find predictions
#   if(N == (object$npar-object$rhoInd)){
#     yinit <- object$Y[(object$nobs - pmax + 1):object$nobs,]
#     yinit <- as.matrix(yinit)
#     muinit <- object$fitted.values[(object$nobs - qmax + 1):object$nobs,]
#     muinit <- as.matrix(muinit)
#
#     newlam <- PoARXPredict(theta, Ynew, Xnew,
#                            object$intercept, object$ylags, object$mulags,
#                            yinit, muinit)
#
#     if(interval & varCond){
#       lamSim <- lapply(seq_len(nrow(coef)), function(j)
#         PoARXPredict(coef[j , ], Ynew, Xnew,
#                      object$intercept, object$ylags, object$mulags,
#                      yinit, muinit))
#       diff <- lapply(lamSim, function(x) x - newlam)
#       diff <- do.call(rbind, diff)
#       quant <- lapply(seq_len(ncol(diff)), function(j){
#         quantile(diff[j,], c(1-alpha, alpha))
#       })
#       quant <- do.call(rbind, quant)
#       lower <- newlam - quant[,1]
#       upper <- newlam - quant[,2]
#       newlam <- rbind(newlam, lower, upper)
#       rownames(newlam) <- c("est", "lower", "upper")
#     }
#   }else{
#     yinit <- object$Y[(object$nobs - pmax + 1):object$nobs,]
#     yinit <- as.matrix(yinit)
#     muinit <- object$fitted.values[(object$nobs - qmax + 1):object$nobs,]
#     muinit <- as.matrix(muinit)
#
#     newlam <- lapply(seq_len(K), function(j){
#       st <- (N*(j-1) + 1)
#       en <- (N*j)
#       PoARXPredict(theta[st:en],
#                    Ynew[ , j , drop = FALSE],
#                    Xnew[ , , j , drop = FALSE],
#                    object$intercept, object$ylags, object$mulags,
#                    yinit[ , j , drop = FALSE],
#                    muinit[ , j , drop = FALSE])
#     })
#     newlam <- do.call(cbind, newlam)
#
#     if(interval & varCond){
#       lamSim <- lapply(seq_len(nrow(coef)), function(i){
#         lam <- lapply(seq_len(K), function(j){
#           st <- (N*(j-1) + 1)
#           en <- (N*j)
#           PoARXPredict(coef[i , st:en],
#                        Ynew[ , j , drop = FALSE],
#                        Xnew[ , , j , drop = FALSE],
#                        object$intercept, object$ylags, object$mulags,
#                        yinit[ , j , drop = FALSE],
#                        muinit[ , j , drop = FALSE])
#         })
#         lam <- do.call(cbind, lam)
#         lam
#       })
#       diff <- lapply(lamSim, function(x) x - newlam)
#       diff <- do.call(rbind, diff)
#       quant <- lapply(seq_len(ncol(diff)), function(j){
#         quantile(diff[j,], c(1-alpha, alpha))
#       })
#       quant <- do.call(rbind, quant)
#       lower <- newlam - quant[,1]
#       upper <- newlam - quant[,2]
#
#       newlam <- lapply(seq_len(K), function(j){
#         cbind(newlam[,j], lower[,j], upper[,j])
#       })
#       newlam <- as.data.frame(do.call(cbind, newlam))
#       names(newlam) <- paste(rep(c("fit", "lwr", "upr"), K),
#                              rep(seq_len(K), times = rep(3,K)),
#                              sep = ".")
#
#     }else{
#       newlam <- as.data.frame(newlam)
#       names(newlam) <- paste("fit", seq_len(K), sep = ".")
#     }
#   }
#   newlam
# }

#' @title Filter mean values for time series following a PoARX model
#'
#' @description
#' Using the observations and the last observed values, use the PoARX filter
#'     equation to predict the unobserved mean values. Where observations
#'     are not found, they are estimated using the relevant mean value.
#'
#' \strong{TODO:} This needs a look. Most of the documentation for this function
#'     seemed to have been copied verbatim from \code{\link{PoARXFilter}},
#'     including the title. If this is intended, we can describe them together.
#'     The descriptions of the parameters below are identical to those in
#'     \code{\link{PoARXFilter}}, except for \code{obs}.
#'
#' @param theta A vector of parameter values.
#' @param horizon A numeric value for the prediction horizon required.
#' @param nts A numeric value detailing the number of time series present.
#' @param obs An optional matrix of discrete-valued time series observations.
#'    Dimensions are as described in \code{\link{PoARXFilter}}. If these are
#'    not provided and more than one-step ahead is required, they will be
#'    estimated using the intensity values.
#' @param exogenous An optional 3-dimensional array containing the values of the
#'    exogenous covariates. Dimensions are as described in
#'    \code{\link{PoARXFilter}}. Must be provided if there are exogenous
#'    covariates in the model.
#' @param intercept A logical indicator for the presence of an intercept in the
#'    model.
#' @param ylags A vector of time lags in the filter equation relating to the
#'    observations.
#' @param mulags A vector of time lags in the filter equation relating to the
#'    intensities.
#' @param yinit A matrix detailing the initialisations of the observations.
#' @param muinit A matrix detailing the initialisations of the means.
#' @param ... currently not used.
#'
#' @return
#'   TODO: Is this exactly as described in \code{\link{PoARXFilter}}?
#'
#' @description
#'   Description is currently identical to \code{\link{PoARXFilter}}
#'
#' @export
PoARXPredict <- function(theta, horizon, nts,
                         obs = NULL, exogenous = NULL,
                         intercept = TRUE, ylags = NULL, mulags = NULL,
                         yinit = NULL, muinit = NULL, ...){
  if(.test_state$get("PoARXFilter") == "0.3.2.9001"){
    ## Perform checks
    if(!is.vector(theta))
      stop("'theta' must be a vector")
    if(!is.array(exogenous))
      stop("'exogenous' must be an array")

    if(!is.null(yinit) & !is.matrix(yinit))
      stop("'yinit' must be a matrix or data.frame")
    if(!is.null(muinit) & !is.matrix(muinit))
      stop("'muinit' must be a matrix or data.frame")

    if(any(length(ylags) == 0))
      ylags <- NULL
    if(any(length(mulags) == 0))
      mulags <- NULL

    if(is.data.frame(yinit))
      yinit <- as.matrix(yinit)
    if(is.data.frame(muinit))
      muinit <- as.matrix(muinit)

    n <- dim(exogenous)[1]
    r <- max(dim(exogenous)[2], 0)
    K <- dim(exogenous)[3]

    w <- as.numeric(intercept)
    P <- length(ylags)
    Q <- length(mulags)
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

    #### Filter the means
    P1 <- min(P,1)
    Q1 <- min(Q,1)
    r1 <- min(r,1)

    Y <- rbind(yinit, obs)
    if(nrow(Y) < n + pmax)
      Y <- rbind(Y, matrix(NA, nrow = n + pmax - nrow(Y), ncol = K))
    lam <- matrix(NA, nrow = qmax + n, ncol = K)
    lam[seq_len(qmax),] <- muinit

    for(t in (seq_len(n) + qmax)){
      newlam <- theta[1]*w +
        (P>0)*colSums(theta[(w+P1):(w+P)] *
                        Y[pmax(t+pmax-qmax-ylags, 0) , , drop = FALSE]) +
        (Q>0)*colSums(theta[(w+P+Q1):(w+P+Q)] *
                        lam[pmax(t-mulags, 0) , , drop = FALSE]) +
        (r>0)*colSums(matrix(theta[(w+P+Q+r1):(w+P+Q+r)] *
                               exogenous[pmax(t-qmax, 0) , , ],
                             ncol = K, byrow = T))
      lam[t,] <- newlam
      obsRow <- Y[t-qmax+pmax,]
      naObsInd <- is.na(obsRow)
      if(any(naObsInd))
        Y[t-qmax+pmax,][naObsInd] <- newlam[naObsInd]
    }
    rownames(lam) <- NULL

    #### Return all observations and mean values
    lam[(qmax + 1):nrow(lam) , , drop = FALSE]
  }else{

    ## - what if there are no exogenous covariates?
    ## - need to replace NAs

    ## To be able to use PoARXPredict, checks must have been performed in
    ##    PoARXFilter - do we need do them again??

#    ## Perform checks
#    if(!is.vector(theta))
#      stop("'theta' must be a vector")
#    if(!is.array(exogenous))
#      stop("'exogenous' must be an array")
#
#    if(!is.null(yinit) & !is.matrix(yinit))
#      stop("'yinit' must be a matrix or data.frame")
#    if(!is.null(muinit) & !is.matrix(muinit))
#      stop("'muinit' must be a matrix or data.frame")

    if(length(ylags) == 0)
      ylags <- NULL
    if(length(mulags) == 0)
      mulags <- NULL

    if(is.data.frame(yinit))
      yinit <- as.matrix(yinit)
    if(is.data.frame(muinit))
      muinit <- as.matrix(muinit)

    w <- as.numeric(intercept)
    P <- length(ylags)
    Q <- length(mulags)
    r <- max(dim(exogenous)[2], 0)

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

    #### Filter the means
    P1 <- min(P,1)
    Q1 <- min(Q,1)
    r1 <- min(r,1)

    Y <- rbind(yinit, obs)
    if(nrow(Y) < horizon + pmax)
      Y <- rbind(Y, matrix(NA, nrow = horizon + pmax - nrow(Y), ncol = nts))

    lam <- matrix(0, nrow = 1, ncol = nts)
    if(!is.null(muinit)){
      lam <- rbind(muinit, lam)
    }

    if(P > 0)
      thetaP <- theta[(w+P1):(w+P)]
    if(Q > 0)
      thetaQ <- theta[(w+P+Q1):(w+P+Q)]
    if(r > 0)
      thetar <- theta[(w+P+Q+r1):(w+P+Q+r)]

    currow <- nrow(lam)
    lam <- rbind(lam, matrix(0, nrow = horizon, ncol = nts))
    pmqy <- pmax - qmax - ylags
    for(t in (seq_len(horizon) + qmax)){
      currow <- currow + 1
      newlam <-
        (if(w > 0) theta[1] else 0) +
        (if(P > 0) .colSums(thetaP * Y[t + pmqy, , drop = FALSE], P, nts)
         else 0) +
        (if(Q > 0) .colSums(thetaQ * lam[t - mulags, , drop = FALSE], Q, nts)
         else 0) +
        (if(r > 0) .colSums(matrix(thetar * exogenous[t - qmax, , ],
                                   ncol = nts, byrow = TRUE), r, nts) else 0)
      lam[currow, ] <- newlam
      naObsInd <- is.na(Y[t-qmax+pmax,])
      if(any(naObsInd))
        Y[t-qmax+pmax,][naObsInd] <- newlam[naObsInd]
    }
    rownames(lam) <- NULL

    #### Return all observations and mean values
    lam[(qmax + 1):nrow(lam) , , drop = FALSE]
  }
}
