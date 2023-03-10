#' @title Marginal PoARX log-likelihood computations
#'
#' @description Marginal log-likelihood computations for PoARX models.
#'
#' @details \code{marginalPoARXLogLikelihood} computes the log-likelihood for
#'     observations from the marginal distributions of PoARX models.
#'
#' @param theta A vector of parameter values.
#' @param obs A matrix of discrete-valued time series observations. Each
#'    column of the matrix represents a separate time series.
#' @param exogenous A 3-dimensional array containing the values of the exogenous
#'    covariates. Each new value in dimension 1 represents another observation
#'    of the time series. Each new vector in dimension 2 corresponds to a
#'    new covariate and each matrix of the third dimension represents a
#'    separate time series.
#' @param intercept A logical indicator for whether an intercept should be
#'    included in the model.
#' @param ylags A vector of time lags in the filter equation relating to the
#'    observations.
#' @param mulags A vector of time lags in the filter equation relating to the
#'    intensities.
#' @param yinit A matrix detailing the initialisations of the observations.
#' @param muinit A matrix detailing the initialisations of the means.
#' @param ... currently not used.
#'
#' @return for \code{marginalPoARXLogLikelihood}, a numeric value representing
#'     the log-likelihood of the observations given the parameters.
#'
#' @export
marginalPoARXLogLikelihood <- function(theta, obs, exogenous = NULL,
                                       intercept = TRUE, ylags = NULL, mulags = NULL,
                                       yinit = NULL, muinit = NULL, ...){
  ## Calculate the means
  obj <- PoARXFilter(theta, obs, exogenous, intercept,
                     ylags, mulags, yinit, muinit)
  ## Use the Poisson model to find the log probability of each observations
  logprobs <- dpois(obs, obj$Lambda$Calculated, log = TRUE)
  ## Due to independence, sum all log probabilities
  sum(logprobs)
}

#' @rdname marginalPoARXLogLikelihood
#'
#' @details \code{marginalPoARXScore} computes marginal PoARX score function,
#'     the score function for observations from the marginal distributions of
#'     PoARX models.
#'
#' @param sameParams A logical argument indicating whether the time series
#'    should be modelled using the same set of PoARX parameter values.
#'
#' @return for \code{marginalPoARXScore}, a numeric vector representing the
#'     score function of the observations given the parameters.
#'
#' @export
marginalPoARXScore <- function(theta, obs, exogenous = NULL,
                               intercept = TRUE, ylags = NULL, mulags = NULL,
                               yinit = NULL, muinit = NULL,
                               sameParams = FALSE, ...){
  score <- .PoARXScore(theta, obs, exogenous,
                       intercept, ylags, mulags,
                       yinit, muinit)
  if(sameParams){
    apply(score, 2, sum)
  }else{
    apply(score, 2:3, sum)
  }
}

#' @title Marginal PoARX score working function
#'
#' @description
#' Find the score function for observations from the marginal distributions of
#'    PoARX models. This returns the score for each individual observation.
#'
#' @inheritParams marginalPoARXLogLikelihood
#'
#' @return An array of dimension \code{(n, N, K)} containing the score vector for
#'    each observation given the parameters.
#'
#' @details \code{n} is the number of observations, \code{N} is the dimension of the
#'    score vector and \code{K} is the number of time series.
#'
#' @keywords internal
#' @rdname dot_PoARXScore
.PoARXScore <- function(theta, obs, exogenous = NULL,
                        intercept = TRUE, ylags = NULL, mulags = NULL,
                        yinit = NULL, muinit = NULL, ...){
    ## Preparation
    K    <- ncol(obs)
    n    <- nrow(obs)
    pmax <- max(ylags, 0)
    qmax <- max(mulags, 0)
    N    <- length(theta)

    ## Calculate the means
    obj <- PoARXFilter(theta, obs, exogenous, intercept,
                       ylags, mulags, yinit, muinit)
    ## Poisson derivative
    A <- (obj$Y$Observed/obj$Lambda$Calculated - 1)
    ## Score
    score <- NULL
    ## Omega derivative
    omega <- if(intercept)
                 matrix(1, ncol = K)
             else NULL

    if(is.null(exogenous)){
        dim1 <- intercept + length(ylags) + length(mulags)
        score <- array(0, dim = c(dim1, K, n))
        for(t in seq_len(n)){
            ## Alpha derivative
            alpha <- obj$Y$All[t+pmax-ylags , , drop = FALSE]
            ## Beta derivative
            beta <- obj$Lambda$All[t+qmax-mulags , , drop = FALSE]
            ## All together
            derivs <- rbind(omega, alpha, beta)

            score[ , , t] <- tcrossprod(derivs, A[t , ])
        }
    }else{
        ## TODO: @jamie: I am not sure about the first dimension of the array,
        ##               so for now I create it inside the loop at the 1st iteration.
        ##               Ugly but should be ok. It would be better to do it here.

        ## dim1 = w + p + q + r
        dim1 <- intercept + length(ylags) + length(mulags) + dim(exogenous)[2]
        score <- array(0, dim = c(dim1, K, n))
        for(t in seq_len(n)){
            ## Alpha derivative
            alpha <- obj$Y$All[t + pmax - ylags , , drop = FALSE]
            ## Beta derivative
            beta <- obj$Lambda$All[t + qmax - mulags , , drop = FALSE]
            ## Eta derivative
            eta <- matrix(exogenous[t , , ], ncol = K)
            ## All together
            derivs <- rbind(omega, alpha, beta, eta)

            if(.test_state$get("PoARXFilter") == "0.3.2"){
                if(t == 1){
                    score <- array(0, dim = c(nrow(derivs), K, n))
                    stopifnot( dim(t(A[t , ] * t(derivs))) == dim(score)[1:2] )
                }
            }
            ## score <- abind(score, t(A[t , ] * t(derivs)), along = 3)
            ##
            ## (AB')' = BA' = tcrossprod(B, A)
            ##     score[ , , t] <- t(A[t , ] * t(derivs))
            score[ , , t] <- tcrossprod(derivs, A[t , ])
        }
    }
    rownames(score) <- NULL
    ## return score output
    aperm(score, c(3,1,2))
}

#' @rdname marginalPoARXLogLikelihood
#'
#' @details \code{marginalPoARXInfoMat} computes the marginal PoARX
#'     information matrix, the information matrix for observations from the
#'     marginal distributions of PoARX models.
#'
#' @param filter A list containing observations and mean values, output from
#'    \code{PoARXFilter}.
#'
#' @return for \code{marginalPoARXInfoMat},
#'    \strong{TODO: Georgi deleted the wrong spec that was here.}
#'
#' @export
marginalPoARXInfoMat <- function(filter, exogenous, intercept = TRUE,
                                 ylags = NULL, mulags = NULL, ...){
    ## Poisson derivative
    A <- (filter$Y$Observed/filter$Lambda$Calculated - 1)
    B <- -filter$Y$Observed/(filter$Lambda$Calculated^2)

    ## Preparation
    K <- ncol(filter$Y$Observed)
    w <- as.numeric(intercept)
    P <- length(ylags)
    Q <- length(mulags)
    r <- max(dim(exogenous)[2], 0)
    N <- w+P+Q+r
    pmax <- max(ylags, 0)
    qmax <- max(mulags, 0)

    ## Omega
    if(intercept){
        omega1 <- matrix(1, ncol = K)
        omega2 <- array(rep(c(0, 1, 0), times = c(1 + P, Q, r)),
          #c(rep(0, 1+P), rep(1, Q), rep(0, r)),
                        dim = c(1, N, K))
    }else{
        omega1 <- NULL
        omega2 <- NULL
    }

    info <- matrix(0, nrow = N, ncol = N)

    if(Q > 0){
        ypast0 <- outer(pmax - mulags, ylags, "-")
        dim(ypast0) <- c(Q, P, K)
        mi0 <- min(ypast0)

        mupast0 <- outer(qmax - mulags, mulags, "-")
        dim(mupast0) <- c(Q, Q, K)

        yarr <- array(0, dim = c(Q, P, K))
        muarr <- array(0, dim = c(Q, Q, K))
    }

    if(is.null(exogenous)){
      if(.test_state$get("PoARXFilter") == "0.3.2"){
        for(t in seq_len(nrow(filter$Y$Observed))){
            ## First derivatives
            alpha1 <- filter$Y$All[t+pmax-ylags , , drop = FALSE]
            beta1 <- filter$Lambda$All[t+qmax-mulags , , drop = FALSE]
            derivs1 <- rbind(omega1, alpha1, beta1)
            derivs1 <- unlist(lapply(seq_len(K), function(j)
                derivs1[ , j] %o% derivs1[ , j]))
            derivs1 <- array(derivs1, dim = c(N, N, K))

            if(Q>0){
                ypast <- pmax(outer(t+pmax-mulags, ylags, "-"), 0)
                ypast <- array(ypast, dim = c(Q, P, K))
                yarr <- array(0, dim = c(Q, P, K))
                yarr[ypast>0] <- filter$Y$All[ypast[ , , 1] , ]

                mupast <- pmax(outer(t+qmax-mulags, mulags, "-"), 0)
                mupast <- array(mupast, dim = c(Q, Q, K))
                muarr <- array(0, dim = c(Q, Q, K))
                muarr[mupast>0] <- filter$Lambda$All[mupast[ , , 1] , ]

                xInd <- pmax(t-mulags, 0)

                alpha2 <- abind(array(0, dim = c(P, w + P, K)),
                                aperm(yarr, c(2,1,3)),
                                array(0, dim = c(P, r, K)),
                                along = 2)
                beta2 <- abind(array(1, dim = c(Q, w, K)),
                               yarr, muarr, along = 2)
                derivs2 <- abind(omega2, alpha2, beta2, along = 1)
            }else{
                derivs2 <- abind(omega2,
                                 array(0, dim = c(P+r, w+P+r, K)),
                                 along = 1)
            }

            info1 <- A[t , ] * aperm(derivs2, c(3,1,2)) +
                B[t , ] * aperm(derivs1, c(3,1,2))
            info <- info + apply(info1, 2:3, sum)
        }
      }else{
        derivs2 <- array(0, dim = c(N, N, K))
        if(intercept){
          derivs2[1 , , ] <- omega2
        }
        yDim <- (w + 1):(w + P)
        lamDim <- (w + P + 1):N # because r = 0
        if(Q > 0){
          ypast0 <- outer(pmax - mulags, ylags, "-")
          dim(ypast0) <- c(Q, P, K)
          mi0 <- min(ypast0)

          mupast0 <- outer(qmax - mulags, mulags, "-")
          dim(mupast0) <- c(Q, Q, K)

          yarr <- array(0, dim = c(Q, P, K))
          muarr <- array(0, dim = c(Q, Q, K))
        }
        for(t in seq_len(nrow(filter$Y$Observed))){
          ## First derivatives
          #alpha1 <- filter$Y$All[t+pmax-ylags , , drop = FALSE]
          #beta1 <- filter$Lambda$All[t+qmax-mulags , , drop = FALSE]
          #derivs1 <- rbind(omega1, alpha1, beta1)
          derivs1 <- rbind(omega1,
                           filter$Y$All[t+pmax-ylags , , drop = FALSE],        # alpha
                           filter$Lambda$All[t+qmax-mulags , , drop = FALSE])  # beta

          ## @jamie @georgi: there must be a function that does this?
          derivs1 <- unlist(lapply(seq_len(K), function(j)
            derivs1[ , j] %o% derivs1[ , j]))
          derivs1 <- array(derivs1, dim = c(N, N, K))

          if(Q>0){
            ##   ypast <- pmax(outer(t + pmax - mulags, ylags, "-"), 0)
            ##   ypast <- array(ypast, dim = c(Q, P, K))
            ypast <- ypast0 <- ypast0 + 1
            ypast[ypast < 0] <- 0
            ## yarr <- array(0, dim = c(Q, P, K))
            yarr[] <- 0
            yarr[ypast > 0] <- filter$Y$All[ypast[ , , 1], ]

            ##   mupast <- pmax(outer(t + qmax - mulags, mulags, "-"), 0)
            ##   mupast <- array(mupast, dim = c(Q, Q, K))
            mupast <- mupast0 <- mupast0 + 1
            mupast[mupast < 0] <- 0
            ## muarr <- array(0, dim = c(Q, Q, K))
            muarr[] <- 0
            muarr[mupast > 0] <- filter$Lambda$All[mupast[ , , 1], ]

            derivs2[lamDim , yDim , ] <-
              derivs2[yDim , lamDim , ] <-
              yarr
            derivs2[lamDim , lamDim , ] <-
              derivs2[lamDim , lamDim , ] <-
              muarr
          }

          info1 <- A[t , ] * aperm(derivs2, c(3,1,2)) +
            B[t , ] * aperm(derivs1, c(3,1,2))
          info <- info + apply(info1, 2:3, sum)
        }
      }
    }else{
      if(.test_state$get("PoARXFilter") == "0.3.2"){
        for(t in seq_len(nrow(filter$Y$Observed))){
            ## First derivatives
            alpha1 <- filter$Y$All[t + pmax - ylags, , drop = FALSE]
            beta1 <- filter$Lambda$All[t + qmax - mulags, , drop = FALSE]
            eta1 <- matrix(exogenous[t , , ], nrow = r)

            derivs1 <- rbind(omega1, alpha1, beta1, eta1)
            derivs1 <- unlist(
                lapply(seq_len(K), function(j) derivs1[ , j] %o% derivs1[ , j]))
            derivs1 <- array(derivs1, dim = c(N, N, K))

            if(Q > 0){
                ##   ypast <- pmax(outer(t + pmax - mulags, ylags, "-"), 0)
                ##   ypast <- array(ypast, dim = c(Q, P, K))
                ypast <- ypast0 <- ypast0 + 1
                ypast[ypast < 0] <- 0
                ## yarr <- array(0, dim = c(Q, P, K))
                yarr[] <- 0 # @georgi: what does this line do?
                yarr[ypast > 0] <- filter$Y$All[ypast[ , , 1], ]

                ##   mupast <- pmax(outer(t + qmax - mulags, mulags, "-"), 0)
                ##   mupast <- array(mupast, dim = c(Q, Q, K))
                mupast <- mupast0 <- mupast0 + 1
                mupast[mupast < 0] <- 0
                ## muarr <- array(0, dim = c(Q, Q, K))
                muarr[] <- 0
                muarr[mupast > 0] <- filter$Lambda$All[mupast[ , , 1], ]

                xInd <- pmax(t - mulags, 0)
                xarr <- array(0, dim = c(Q, r, K))
                xarr[xInd > 0, , ] <- exogenous[xInd, , ]

                ## TODO: (@jamie:) this probably can be sped up substantially
                ##       but I need the maths.
                alpha2 <- abind(array(0, dim = c(P, w + P, K)),
                                aperm(yarr, c(2,1,3)),
                                array(0, dim = c(P, r, K)),
                                along = 2)
                beta2 <- abind(array(1, dim = c(Q, w, K)),
                               yarr, muarr, xarr,
                               along = 2)
                eta2 <- abind(array(0, dim = c(r, w + P, K)),
                              aperm(xarr, c(2,1,3)),
                              array(0, dim = c(r, r, K)),
                              along = 2)
                derivs2 <- abind(omega2, alpha2, beta2, eta2, along = 1)
            }else{
                derivs2 <- abind(omega2, array(0, dim = c(P+r, w+P+r, K)),
                                 along = 1)
            }

            info1 <- A[t , ] * aperm(derivs2, c(3,1,2)) +
                     B[t , ] * aperm(derivs1, c(3,1,2))
            info <- info + apply(info1, 2:3, sum)
        }
      }else{
        derivs2 <- array(0, dim = c(N, N, K))
        if(intercept){
          derivs2[ 1 , , ] <- derivs2[ , 1 , ] <- omega2
        }
        yDim <- (w + 1):(w + P)
        lamDim <- (w + P + 1):(w + P + Q)
        exDim <- (w + P + Q + 1):N
        for(t in seq_len(nrow(filter$Y$Observed))){
          ## First derivatives
          #alpha1 <- filter$Y$All[t + pmax - ylags, , drop = FALSE]
          #beta1 <- filter$Lambda$All[t + qmax - mulags, , drop = FALSE]
          #eta1 <- matrix(exogenous[t , , ], nrow = r)
          #derivs1 <- rbind(omega1, alpha1, beta1, eta1)

          derivs1 <- rbind(omega1,
                           filter$Y$All[t+pmax-ylags , , drop = FALSE],        # alpha
                           filter$Lambda$All[t+qmax-mulags , , drop = FALSE],  # beta
                           matrix(exogenous[t , , ], nrow = r))                # eta

          ## @jamie @georgi: must be a better way to do this?
          derivs1 <- unlist(
            lapply(seq_len(K), function(j) derivs1[ , j] %o% derivs1[ , j]))
          derivs1 <- array(derivs1, dim = c(N, N, K))

          if(Q > 0){
            ##   ypast <- pmax(outer(t + pmax - mulags, ylags, "-"), 0)
            ##   ypast <- array(ypast, dim = c(Q, P, K))
            ypast <- ypast0 <- ypast0 + 1
            ypast[ypast < 0] <- 0
            ## yarr <- array(0, dim = c(Q, P, K))
            yarr[] <- 0 # @georgi: what does this line do?
            yarr[ypast > 0] <- filter$Y$All[ypast[ , , 1], ]

            ##   mupast <- pmax(outer(t + qmax - mulags, mulags, "-"), 0)
            ##   mupast <- array(mupast, dim = c(Q, Q, K))
            mupast <- mupast0 <- mupast0 + 1
            mupast[mupast < 0] <- 0
            ## muarr <- array(0, dim = c(Q, Q, K))
            muarr[] <- 0
            muarr[mupast > 0] <- filter$Lambda$All[mupast[ , , 1], ]

            xInd <- pmax(t - mulags, 0)
            xarr <- array(0, dim = c(Q, r, K))
            xarr[xInd > 0, , ] <- exogenous[xInd, , ]

            derivs2[lamDim , yDim , ] <-
              derivs2[yDim , lamDim , ] <-
              yarr
            derivs2[lamDim , lamDim , ] <-
              derivs2[lamDim , lamDim , ] <-
              muarr
            derivs2[lamDim , exDim , ] <-
              derivs2[exDim , lamDim , ] <-
              xarr
          }
          ## If there's no past mean values, the second deriv matrix is 0

          info1 <- A[t , ] * aperm(derivs2, c(3,1,2)) +
            B[t , ] * aperm(derivs1, c(3,1,2))
          info <- info + apply(info1, 2:3, sum)
        }
      }
    }
    - info
}

#' @title Full PoARX log-likelihood
#'
#' @description
#' Find the log-likelihood for observations from a PoARX model.
#'
#' @param vartheta A vector of parameter values, including all marginal parameters
#'    and any dependence parameters.
#' @param obs A matrix of discrete-valued time series observations. Each
#'    column of the matrix represents a separate time series.
#' @param xlist A list of matrices containing the values of the exogenous
#'    covariates. Each new value in dimension 1 represents another observation
#'    of the time series. Each new vector in dimension 2 corresponds to a
#'    separate new covariate. Each component of the list represents a
#'    time series.
#' @param modelSummary A data.frame summarising the features of the PoARX models.
#' @param ylagList A list of vectors of time lags in the filter equation
#'    relating to the observations.
#' @param mulagList A list of vectors of time lags in the filter equation
#'    relating to the intensities.
#' @param yinitList A list of matrices detailing the initialisations of the
#'    observations.
#' @param muinitList A list of matrices detailing the initialisations of the means.
#' @param indep A logical argument indicating whether the PoARX models are
#'    coupled using Frank's copula or the independence copula.
#' @param sameParams A logical argument indicating whether the time series
#'    should be modelled using the same set of PoARX parameter values.
#' @param zeroFix A numeric value to transform any 0-probability values to.
#' @param ... currently not used.
#'
#' @return A numeric value representing the log-likelihood of the observations
#'    given the parameters.
#'
#' @export
fullPoARXLogLikelihood <- function(vartheta, obs, xlist,
                                   modelSummary, ylagList, mulagList,
                                   yinitList, muinitList,
                                   indep = FALSE, sameParams = FALSE,
                                   zeroFix = 1e-35, ...){
    ## Prelims
    K <- ncol(obs)
    Kseq <- seq_len(K)
    N <- rowSums(modelSummary)

    ## Sort marginal and dependence parameters out
    if(indep){
        theta <- vartheta
    }else{
        rho <- vartheta[length(vartheta)]
        theta <- vartheta[-length(vartheta)]
    }
    ## Calculate the means
    intercept <- as.logical(modelSummary$w)
    r <- modelSummary$r

    if(sameParams){
        Xarr <- if(r > 0)
                    array(unlist(xlist), dim = c(nrow(obs), r, K))
                else NULL

        yinit <- do.call(cbind, yinitList)
        if(nrow(yinit) == 0)
            yinit <- NULL
        muinit <- do.call(cbind, muinitList)
        if(nrow(muinit) == 0)
            muinit <- NULL

        obj <- PoARXFilter(theta, obs, Xarr, intercept,
                           ylagList[[1]], mulagList[[1]], yinit, muinit
                           )$Lambda$Calculated
    }else{
        ## parameters
        mark <- rep(Kseq, times = N)
        thetaj <- split(theta, mark)
        names(thetaj) <- NULL

        fun <- function(j){
            ob <- as.matrix(obs[ , j, drop = FALSE])
            ## modelSummary[j, "r"]
            ex <- if(r[j] > 0)
                      array(xlist[[j]], dim = c(dim(xlist[[j]]),1))
                  else NULL
            obsinit <- yinitList[[j]]
            if(nrow(obsinit) == 0)
                obsinit <- NULL
            meaninit <- muinitList[[j]]
            if(nrow(meaninit) == 0)
                meaninit <- NULL

            PoARXFilter(theta = thetaj[[j]], ob, ex, intercept[j],
                        ylagList[[j]], mulagList[[j]], obsinit, meaninit
                        )$Lambda$Calculated
        }
        obj <- lapply(seq_len(K), fun)
        obj <- do.call(cbind, obj)
    }

    ## Find log-probabilities
    logprobs <- if(indep)
                    dpois(obs, obj, log = TRUE)
                else
                    log(dFranksCopulaPoisson(obs, obj, rho))

    ## Fit lower bound
    logprobs <- pmax(logprobs, log(zeroFix))

    ## Sum all log probabilities
    sum(logprobs)
}
