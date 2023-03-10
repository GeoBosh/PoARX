#' @title Fit a PoARX model
#'
#' @description Fit a PoARX models to \code{K} time series using the method of
#'     inference functions \insertCite{Joe2005}{PoARX}.
#'
#' @param formula A formula object representing the response values and exogenous
#'    covariates. Different time series can be separated using \code{|}.
#' @param data,subset arguments controlling formula processing via \code{model.frame}.
#' @param init An initial vector of parameter values.
#' @param ylags A vector or matrix of time lags in the filter equation relating to the
#'    observations. If a matrix, the number of columns should correspond to the number
#'    of time series and each column gives the lags of the corresponding time series.
#'    When the lengths of lags do not match, use zeroes to fill the matrix.
#' @param mulags A vector or matrix of time lags in the filter equation relating to the
#'    intensities. If a matrix, the number of columns should correspond to the number
#'    of time series and each column gives the lags of the corresponding time series.
#'    When the lengths of lags do not match, use zeroes to fill the matrix.
#' @param yinit A matrix detailing the initialisations of the observations. If the
#'    maximum lags do not match, fill the required entries with zeroes.
#' @param muinit A matrix detailing the initialisations of the means. If the
#'    maximum lags do not match, fill the required entries with zeroes.
#' @param indep A logical argument indicating whether the PoARX models are
#'    coupled using Frank's copula or the independence copula.
#' @param sameParams A logical argument indicating whether the time series
#'    should be modelled using the same set of PoARX parameter values.
#' @param maxRho A value indicating the largest absolute value that the dependence
#'    parameter of Frank's copula can take. For use in \code{optimise}.
#' @param zeroFix A numeric value to transform any 0-probability values to.
#' @param ... currently not used.
#'
#' @return An object of class 'PoARX'. The structure of the object is internal
#'     and subject to change. Currently it contains the following components:
#'
#'     \item{formula}{model formula, an object from class \code{"Formula".}}
#'     \item{data}{the model frame (?)}
#'     \item{X}{}
#'     \item{Y}{}
#'     \item{nobs, nts}{The number of observations and time series, respectively.}
#'     \item{npar}{}
#'     \item{modelSummary}{}
#'     \item{ylags}{}
#'     \item{ylagList}{}
#'     \item{mulags}{}
#'     \item{mulagList}{}
#'     \item{sameParams}{}
#'     \item{rhoInd}{Information regarding the format of the PoARX model.}
#'     \item{yinit, muinit}{Information regarding the format of the PoARX model.}
#'     \item{rhoValue}{}
#'     \item{theta}{The estimated parameters.}
#'     \item{logLik}{The log-likelihood value for the parameters given the data.}
#'     \item{varCovar}{The variance-covariance matrix of the parameters.}
#'     \item{fitted.values}{}
#'     \item{fitted.probs}{}
#'     \item{intercept}{\strong{TODO: is this really present?}}
#'
#' @details
#' A time series, \eqn{Y_t}, following a PoARX model obeys the following equations:
#'       \deqn{Y_t = N_t(\lambda_t),}
#'       \deqn{\lambda_t = \omega + \sum \alpha_l Y_{t-j_l} + \sum \beta_l
#'                \lambda_{t-j_l} + \eta x_{t-1},}
#'    where \eqn{N_t(.)} represents a Poisson process of unit intensity, so that
#'    \eqn{N_t(\lambda_t)} represents a Poisson process of intensity \eqn{\lambda_t}.
#'    \eqn{Y_t} represents the observations and \eqn{x_t} the covariates at time \eqn{t}.
#'    \eqn{\omega}, {\eqn{\alpha_i}}, {\eqn{\beta_i}} and \eqn{\eta} are all positive
#'    parameters that must be supplied. An additional constraint on the parameters is:
#'       \deqn{\sum (\alpha_i + \beta_i) < 1,}
#'    which means that the model satisfies the conditions for stationarity and
#'    ergodicity. Multiple time series are modelled independently by this fit
#'    procedure. The number of time series is denoted by \code{K}, the number of past
#'    observations by \code{p}, the number of past means by \code{q}, the number of
#'    exogenous covariates by \code{r}, and \code{w} is an indicator of an intercept
#'    in the model.
#'
#' @references
#'    \insertAllCited{}
#'
#' @examples
#' #### examples here!
#'
#' @export
fitPoARX <- function(formula, data, subset, init,
                     ylags = NULL, mulags = NULL,
                     yinit = NULL, muinit = NULL,
                     indep = FALSE, sameParams = FALSE,
                     maxRho = 1000, zeroFix = 1e-50, ...){
    ## formula
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset"), names(mf), 0)
    mfA <- mf[c(1, m)]
    fo <- Formula(formula)
    mfA[[1]] <- as.name("model.frame")
    mfA$formula <- fo
    mfA <- eval(mfA, parent.frame())

    ## preliminary definitions
    n <- nrow(mfA)
    K <- length(fo)[1]
    if(K == 1){
        indep = TRUE
        sameParams = TRUE
    }
    Kseq <- seq_len(K)

    #### data, number of exogenous covariates and intercept terms
    Y <- model.part(fo, data = mfA, lhs = Kseq)
    ynames <- names(Y)
    xlist <- lapply(seq_len(K),
                    function(i) model.matrix(fo, data = mfA, rhs = i))
    xnames <- sapply(xlist, function(z) names(z))
    w <- lapply(xlist, function(Xi) "(Intercept)" %in% colnames(Xi))
    w <- as.numeric(unlist(w))
    xlist <- lapply(xlist, function(Xi)
                               Xi[ , colnames(Xi) != "(Intercept)", drop = FALSE])
    r <- unlist(lapply(xlist, ncol))
    X <- do.call(cbind, xlist)

    ## number of lags, exogenous covariates, npar
    if(is.vector(ylags))
        ylags <- matrix(rep(ylags, K), ncol = K)

    if(is.null(ylags)){
        ylagList <- NULL
        P <- 0
    }else{
        ylagList <- split(t(ylags), Kseq)
        ylagList <- lapply(ylagList, function(x) x[x > 0])
        names(ylagList) <- NULL
        P <- sapply(ylagList, function(x) sum(x > 0))
    }
    if(is.vector(mulags))
        mulags <- matrix(rep(mulags, K), ncol = K)

    if(is.null(mulags)){
        mulagList <- NULL
        Q <- 0
    }else{
        mulagList <- split(t(mulags), Kseq)
        mulagList <- lapply(mulagList, function(x) x[x > 0])
        names(mulagList) <- NULL
        Q <- sapply(mulagList, function(x) sum(x > 0))
    }

    objY <- as.data.frame(Y)

    ## Model summary
    modelSummary <- data.frame(w = w, P = P, Q = Q, r = r)
    ndup <- sum(duplicated(modelSummary))
    if(sameParams){
        if(ndup != (K-1))
            stop("If models are to have same parameters, they must have the same make-up!")
        else
            modelSummary <- modelSummary[!duplicated(modelSummary), , drop = FALSE]
    }
    M <- nrow(modelSummary)
    Mseq <- seq_len(M)
    N <- rowSums(modelSummary)
    .parNamesVecBuild <- function(rowNo, modelSummary){
        w <- modelSummary[rowNo, "w"]
        P <- modelSummary[rowNo, "P"]
        Q <- modelSummary[rowNo, "Q"]
        r <- modelSummary[rowNo, "r"]
        parNames <- NULL
        if(w > 0)
            parNames <- c(parNames, "omega")
        if(P > 0)
            parNames <- c(parNames, paste0("alpha", seq_len(P)))
        if(Q > 0)
            parNames <- c(parNames, paste0("beta", seq_len(Q)))
        if(r > 0)
            parNames <- c(parNames, paste0("eta", seq_len(r)))
        parNames
    }
    parNames <- lapply(Mseq, .parNamesVecBuild, modelSummary = modelSummary)
    if(M>1)
        parNames <- lapply(Mseq, function(i) paste0(parNames[[i]], ".", i))
    parNames <- unlist(parNames)

    nReq <- sum(N) + !indep
    if(nReq == 0)
        stop("No model to fit")

    if(length(init) != nReq)
        stop(paste0("length of 'init', (", length(init),
                    ") does not match the number of required parameters, (",
                    nReq, ")"))

    ## initialisations
    if(is.null(yinit) && !is.null(ylags)){
        yinit <- data.frame(matrix(0, max(ylags), K))
        names(yinit) <- ynames
    }
    if(is.matrix(yinit)){
        yinit <- data.frame(yinit)
        names(yinit) <- ynames
    }
    if(is.null(muinit) && !is.null(mulags)){
        muinit <- data.frame(matrix(0, max(mulags), K))
        names(muinit) <- paste0("lam", seq_len(K))
    }
    if(is.matrix(muinit)){
        muinit <- data.frame(muinit)
        names(muinit) <- paste0("lam", seq_len(K))
    }

    ## exogenous covariates
    if(length(xlist) == 1 && K > 1)
        xlist <- lapply(seq_len(K), function(j) xlist[[1]])

    if(length(xlist) != K)
        stop("lengths of lhs and rhs of formula do not match")

    ## initital value lists
    yinitList <- NULL
    if(!is.null(yinit)){
        .yinitTrim <- function(j){
            pmax1 <- max(ylagList[[j]], 0)
            obsinit <- yinit[ , j, drop = FALSE]
            n1 <- nrow(obsinit)
            if(n1 > pmax1){
                desRows <- (n1):(n1-pmax1+1)
                obsinit <- obsinit[desRows, , drop = FALSE]
            }
            if(pmax1 == 0)
                obsinit <- matrix(0, nrow = 0, ncol = 1)
            obsinit
        }
        yinitList <- lapply(Kseq, .yinitTrim)
        names(yinitList) <- NULL
    }

    muinitList <- NULL
    if(!is.null(muinit)){
        .muinitTrim <- function(j){
            qmax1 <- max(mulagList[[j]], 0)
            meaninit <- muinit[ , j, drop = FALSE]
            n1 <- nrow(meaninit)
            if(n1 > qmax1 && qmax1 > 0){
                desRows <- (n1):(n1-qmax1+1)
                meaninit <- meaninit[desRows, , drop = FALSE]
            }
            if(qmax1 == 0){
                meaninit <- matrix(0, nrow = 0, ncol = 1)
            }
            meaninit
        }
        muinitList <- lapply(Kseq, .muinitTrim)
        names(muinitList) <- NULL
    }

    if(sameParams){
        ## model make up
        intercept <- as.logical(modelSummary$w)
        ylags1 <- ylags[, 1]
        mulags1 <- mulags[, 1]
        r <- modelSummary$r

        ## exogenos array
        Xarr <- if(r > 0)
                    Xarr <- array(unlist(xlist), dim = c(n, r, K))
                else NULL

        ## constraints
        constraints <- PoARXConstraints(M = 1, r = r, intercept = intercept,
                                        ylags = ylags1, mulags = mulags1)

        ## Ensure y, yinit, and muinit are matrices
        Ymat <- as.matrix(Y)

        yinmat <- muinmat <- NULL
        if(!is.null(yinit))
            yinmat <- as.matrix(yinit)
        if(!is.null(muinit))
            muinmat <- as.matrix(muinit)

        ## Define marginal parameters
        theta <- if(indep)
                     init
                 else
                     ## assuming there's one dependence parameter
                     theta <- init[-length(init)]

        ## fit marginal parameters
        fu <- function(j){
            fit <- constrOptim(
                theta      = theta,
                f          = marginalPoARXLogLikelihood,
                grad       = marginalPoARXScore,
                ui         = constraints$ui,
                ci         = constraints$ci,
                obs        = Ymat,
                exogenous  = Xarr,
                intercept  = intercept,
                ylags      = ylags1,
                mulags     = mulags1,
                yinit      = yinmat,
                muinit     = muinmat,
                sameParams = sameParams,
                control    = list(fnscale = -1))

            filter <- PoARXFilter(fit$par, Ymat, Xarr, intercept,
                                  ylags1, mulags1, yinmat, muinmat)

            score <- .PoARXScore(fit$par, Ymat, Xarr, intercept,
                                 ylags1, mulags1, yinmat, muinmat, sameParams)

            infomat <- marginalPoARXInfoMat(filter, Xarr, intercept,
                                            ylags1, mulags1)

            list(Lambda = filter$Lambda$Calculated, Theta  = fit$par,
                 Value  = fit$value, Score  = score, InfoMat = infomat )
        }

        margFit <- list(fu(1))
    }else{
        ## We will fit each time series separately in lists
        ## r is a vector length K
        r <- modelSummary$r

        ## ylist
        ylist <- lapply(seq_len(K),
                        function(j) model.part(fo, data = mfA, lhs = j))

        ## initial theta sort
        theta0 <- init
        if(!indep)
            theta0 <- theta0[-length(theta0)]
        mark <- rep(seq_len(K), times = N)
        thetaj <- split(theta0, mark)
        names(thetaj) <- NULL

        ## fit marginal parameters
        fu <- function(j){
            ob <- as.matrix(ylist[[j]])
            intercept <- as.logical(w[j])

            ylags <- ylagList[[j]]
            if(length(ylags) == 0){
                ylags <- NULL
                obsinit <- NULL
            }else{
                obsinit <- as.matrix(yinitList[[j]])
                ## if(any(nrow(obsinit) == 0)) # Why any here? isn't nrow(obsinit) a scalar?
                if(length(obsinit) == 0)
                    obsinit <- NULL
            }

            mulags <- mulagList[[j]]
            if(length(mulags) == 0){
                mulags <- NULL
                meaninit <- NULL
            }else{
                meaninit <- as.matrix(muinitList[[j]])
                ##if(any(nrow(meaninit) == 0)) # Why any here? isn't nrow(meaninit) a scalar?
                if(length(meaninit) == 0)
                    meaninit <- NULL
            }

            ex <- if(r[j] > 0)
                      array(xlist[[j]], dim = c(dim(xlist[[j]]),1))
                  else NULL

            constraints <- PoARXConstraints(M = 1, r = r[j],
                                            intercept = intercept,
                                            ylags = ylags, mulags = mulags)
            fit <- constrOptim(
                theta     = thetaj[[j]],
                f         = marginalPoARXLogLikelihood,
                grad      = marginalPoARXScore,
                ui        = constraints$ui,
                ci        = constraints$ci,
                obs       = ob,
                exogenous = ex,
                intercept = intercept,
                ylags     = ylags,
                mulags    = mulags,
                yinit     = obsinit,
                muinit    = meaninit,
                control   = list(fnscale = -1))

            filter <- PoARXFilter(fit$par, ob, ex, intercept,
                                  ylags, mulags, obsinit, meaninit)
            score <- .PoARXScore(fit$par, ob, ex, intercept,
                                 ylags, mulags, obsinit, meaninit)
            infomat <- marginalPoARXInfoMat(filter, ex, intercept, ylags, mulags)

            list(Lambda = filter$Lambda$Calculated, Theta  = fit$par,
                 Value  = fit$value, Score  = score, InfoMat = infomat )
        }
        margFit <- lapply(seq_len(K), fu)
    }
    ## Collect mean values
    means <- lapply(Mseq, function(j) margFit[[j]]$Lambda)
    means <- do.call(cbind, means)
    colnames(means) <- names(muinit)

    ## Collect marginal parameter values
    params <- lapply(Mseq, function(j) margFit[[j]]$Theta)
    params <- do.call(c, params)
    names(params) <- parNames

    ## Add rho to parameter names
    if(!indep)
        parNames <- c(parNames, "rho")

    if(indep){
        ## V-matrix is just inverse of info matrices in independent case
        jthInfoMat <- function(j){
            mat <- tryCatch(
                solve(margFit[[j]]$InfoMat),
                error = function(e){
                    if(sameParams){
                        Ymat <- as.matrix(Y)
                        yinmat <- if(is.null(yinit)) NULL else as.matrix(yinit)
                        muinmat <- if(is.null(muinit)) NULL else as.matrix(muinit)
                        hess1 <- - hessian(marginalPoARXLogLikelihood,
                                           x = margFit[[j]]$Theta,
                                           obs = Ymat, exogenous = Xarr,
                                           intercept = intercept,
                                           ylags = ylags, mulags = mulags,
                                           yinit = yinmat, muinit = muinmat)
                    }else{
                        ob <- as.matrix(ylist[[j]])

                        ex <- if(r[j] > 0)
                                  array(xlist[[j]], dim = c(dim(xlist[[j]]),1))
                              else NULL

                        if(length(yinitList[[j]]) > 0){
                          obsinit <- as.matrix(yinitList[[j]])
                        }else{
                          obsinit <- NULL
                        }

                        if(length(muinitList[[j]]) > 0){
                          meaninit <- as.matrix(muinitList[[j]])
                        }else{
                          meaninit <- NULL
                        }

                        hess1 <- - hessian(marginalPoARXLogLikelihood,
                                           x = margFit[[j]]$Theta,
                                           obs = ob, exogenous = ex,
                                           intercept = intercept,
                                           ylags = ylags, mulags = mulags,
                                           yinit = obsinit, muinit = meaninit)
                    }
                    tryCatch(solve(hess1), function(e) NULL)
                })
            if(!is.null(mat))
                as.matrix(nearPD(mat)$mat)
            else
                NULL
        }

        Vmat <- lapply(Mseq, jthInfoMat)
        nullInd <- sapply(Vmat, is.null)
        Vmat <- if(any(nullInd)){
                    message("Variance matrix cannot be calculated")
                    NULL
                }else
                    as.matrix(do.call(bdiag, Vmat))

        probs <- dpois(as.matrix(objY), means)
    }else{
        ## fit rho
        obs <- as.matrix(objY)
        if(K==2){
            rho_lower <- - maxRho
            rho_hessian <- FALSE
        }else{
            rho_lower <- 0
            rho_hessian <- TRUE
        }

        rhoFit <- optim(
            par     = init[length(init)],
            fn      = rhoFCLogLikelihood,
            gr      = rhoFCScore,
            method  = "Brent",
            lower   = rho_lower,
            upper   = maxRho,
            obs     = obs,
            means   = means,
            control = list(fnscale = -1),
            hessian = rho_hessian )
        params <- c(params, rho = rhoFit$par)

        ## D-matrix
        Idiag <- lapply(Mseq, function(j) margFit[[j]]$InfoMat )
        ## check all can be inverted
        cond <- sapply(Mseq,
                       function(j){
                           tryCatch({
                               inv <- solve(Idiag[[j]])
                               TRUE
                           }, error = function(e) FALSE)
                       })
        ## if not, estimate numerically
        if(!all(cond)){
            jth_hess <- function(j){
                ob <- as.matrix(ylist[[j]])
                ex <- if(r > 0)  # r[j] ???  !!!
                          array(xlist[[j]], dim = c(dim(xlist[[j]]),1))
                      else NULL


                if(length(yinitList[[j]]) > 0){
                  obsinit <- as.matrix(yinitList[[j]])
                }else{
                  obsinit <- NULL
                }

                if(length(muinitList[[j]]) > 0){
                  meaninit <- as.matrix(muinitList[[j]])
                }else{
                  meaninit <- NULL
                }

                - hessian(marginalPoARXLogLikelihood,
                          x = margFit[[j]]$Theta,
                          obs = ob, exogenous = ex,
                          intercept = intercept,
                          ylags = ylags, mulags = mulags,
                          yinit = obsinit, muinit = meaninit)
            }

            Idiag[!cond] <- lapply(Mseq[!cond], jth_hess)
        }
        Idiag <- as.matrix(do.call(bdiag, Idiag))

        .findBottomRow <- function(zeroFix){
            Irhoj <- .rhothetaFCInfoMat(params, obs, means, xlist,
                                        modelSummary, ylagList, mulagList,
                                        yinitList, muinitList,
                                        sameParams, zeroFix)
            Irhorho <- .rhoFCInfoMat(rhoFit$par, obs = obs, means = means,
                                     zeroFix = zeroFix)
            c(Irhoj, Irhorho)
        }

        matUndefined <- TRUE
        Dinv <- NULL
        nsteps <- 0
        while(matUndefined && nsteps < 10){
            nsteps <- nsteps + 1
            bottomRow <- .findBottomRow(zeroFix)
            minusDmat <- rbind(cbind(Idiag, 0), bottomRow)
            dim1 <- ncol(minusDmat)
            matUndefined <- FALSE
            Dinv <- tryCatch(solve(minusDmat),
                             error = function(e){
                                 message("'zeroFix' too small, increasing by a factor of 1e5")
                                 zeroFix <<- zeroFix * 1e5
                                 NULL
                             })
            matUndefined <- is.null(Dinv) | any(is.na(Dinv))
        }

        if(!matUndefined){
            ## M-matrix
            g <- lapply(Mseq, function(j) margFit[[j]]$Score[ , , 1] )
            g <- do.call(cbind, g)
            gp <- rhoFCScore(rhoFit$par, obs, means, zeroFix)
            Mmat <- as.matrix(bdiag(cov(g), var(gp)))

            ## V-matrix
            Vmat <- Dinv * Mmat * Dinv
            Vmat <- as.matrix(nearPD(Vmat)$mat)
        }else{
            message("Variance matrix cannot be calculated")
            Vmat <- NULL
        }
        probs <- dFranksCopulaPoisson(objY, means, rhoFit$par)
    }
    ## Vmat is the variance matrix for sqrt(n) * vartheta, so divide by n
    if(!is.null(Vmat)){
        Vmat <- Vmat/n
        rownames(Vmat) <- colnames(Vmat) <- parNames
    }
    names(probs) <- NULL
    probs <- matrix(probs, nrow = n)

    probs_df <- as.data.frame(probs)
    names(probs_df) <- if(ncol(probs_df) > 1)
                           paste0("probs", seq_len(K))
                       else "probs"

    if(indep){
      structure(list(
        formula       = fo,
        data          = mfA,
        Y             = objY,
        X             = as.data.frame(X),
        nobs          = n,
        nts           = K,
        npar          = sum(N) + (!indep),
        modelSummary  = modelSummary,
        ylags         = ylags,
        ylagList      = ylagList,
        mulags        = mulags,
        mulagList     = mulagList,
        sameParams    = sameParams,
        rhoInd        = !indep,
        yinit         = yinit,
        muinit        = muinit,
        rhoValue      = NULL,
        theta         = params,
        logLik        = sum(pmax(log(probs), log(zeroFix))),
        varCovar      = Vmat,
        fitted.values = as.data.frame(means),
        fitted.probs  = probs_df
      ), class = "PoARX")
    }else{
      structure(list(
        formula       = fo,
        data          = mfA,
        Y             = objY,
        X             = as.data.frame(X),
        nobs          = n,
        nts           = K,
        npar          = sum(N) + (!indep),
        modelSummary  = modelSummary,
        ylags         = ylags,
        ylagList      = ylagList,
        mulags        = mulags,
        mulagList     = mulagList,
        sameParams    = sameParams,
        rhoInd        = !indep,
        yinit         = yinit,
        muinit        = muinit,
        rhoValue      = rhoFit$par,
        theta         = params,
        logLik        = sum(pmax(log(probs), log(zeroFix))),
        varCovar      = Vmat,
        fitted.values = as.data.frame(means),
        fitted.probs  = probs_df
      ), class = "PoARX")
    }
}
