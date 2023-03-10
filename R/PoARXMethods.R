#' @title Methods for PoARX objects
#'
#' @description
#' Standard methods for PoARX objects.
#'
#' @param object,x A \code{"PoARX"} object from \code{PoARXFit}.
#' @param ... Further arguments.
#'
#' @details
#' Objects from class \code{"PoARX"} represent a fitted PoARX model, independent or
#'    dependent, as fitted by \code{indPoARXFit} or \code{depPoARXFit}, respectively.
#'
#' @examples
#' #### Examples go here!
#'
#' @keywords methods
#' @name PoARX_methods
NULL

#' @rdname PoARX_methods
#'
#' @export
AIC.PoARX <- function(object, ...){
  2*(object$npar - object$logLik)
}

#' @rdname PoARX_methods
#'
#' @export
coef.PoARX <- function(object, ...){
  object$theta
}

#' @rdname PoARX_methods
#'
#' @export
df.residual.PoARX <- function(object, ...){
  if(object$sameParams){
    if(object$rhoInd){
      message("rho has different number of degrees of freedom to other parameters")
      c(object$nobs*object$nts - object$npar,
        object$nobs - object$npar)
    }else{
      object$nobs*object$nts - object$npar
    }
  }else{
    object$nobs-object$npar
  }
}

#' @rdname PoARX_methods
#'
#' @export
fitted.PoARX <- function(object, ...){
  object$fitted.values
}

#' @rdname PoARX_methods
#'
#' @export
logLik.PoARX <- function(object, ...){
  object$logLik
}

#' @rdname PoARX_methods
#'
#' @export
nobs.PoARX <- function(object, ...){
  object$nobs
}

#' @rdname PoARX_methods
#'
#' @export
print.PoARX <- function(x, ...){
  w <- x$modelSummary$w
  nMod <- length(w)
  ylags <- x$ylags
  if(is.null(ylags)){
    ylags <- matrix(0, ncol = nMod)
  }else if(is.vector(ylags)){
    ylags <- matrix(ylags, ncol = 1)
  }
  mulags <- x$mulags
  if(is.null(mulags)){
    mulags <- matrix(0, ncol = nMod)
  }else if(is.vector(mulags)){
    mulags <- matrix(mulags, ncol = 1)
  }
  r <- x$modelSummary$r

  for(i in seq_len(nrow(x$modelSummary))){
    if(nrow(x$modelSummary) > 1){
      cat(paste0("Model ", i, ":\n"))
    }else{
      cat("Model:\n")
    }
    cat(paste("intercept:            ", w[i], "\n"))
    cat(paste("observation lags:     ", paste(ylags[,i], collapse = ", "), "\n"))
    cat(paste("mean lags:            ", paste(mulags[,i], collapse = ", "), "\n"))
    cat(paste("exogenous covariates: ", r[i], "\n"))
    cat("\n")
  }
  cat("Coefficients:\n")
  print.default(format(x$theta, digits = 5),
                print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @rdname PoARX_methods
#'
#' @export
residuals.PoARX <- function(object, ...){
  object$Y - object$fitted.values
}

#' @rdname PoARX_methods
#'
#' @export
se.coef.PoARX <- function(object, ...){
  if(is.null(object$varCovar)){
    message("Standard errors not available")
    NULL
  }else{
    sqrt(diag(object$varCovar))
  }
}

#' @rdname PoARX_methods
#'
#' @export
summary.PoARX <- function(object, ...){
  ## coefficient matrix
  if(object$sameParams){
    if(object$rhoInd){
      n1 <- object$nobs*object$nts - object$npar
      n2 <- object$nobs - object$npar
    }else{
      n1 <- object$nobs*object$nts - object$npar
    }
  }else{
    n1 <- object$nobs-object$npar
  }
  npar <- length(object$theta)
  if(is.null(object$varCovar)){
    s1 <- t1 <- p1 <- stars <- rep(NA, npar)
  }else{
    s1 <- signif(sqrt(diag(object$varCovar)), 5)
    t1 <- object$theta/sqrt(diag(object$varCovar))
    p1 <- 2*(1-pt(abs(t1), n1))
    if(exists("n2")){
      p1[npar] <- 2*(1-pt(abs(t1[npar]), n2))
    }
    p1 <- signif(p1, 5)
    stars <- rep(" ", npar)
    stars <- ifelse(p1 < 0.1, ".", stars)
    stars <- ifelse(p1 < 0.05, "*", stars)
    stars <- ifelse(p1 < 0.01, "**", stars)
    stars <- ifelse(p1 < 0.001, "***", stars)
    p1[p1 == 0] <- "<= 1e-16"
  }
  obj <- matrix(c(names(object$theta),
                  signif(object$theta, 5),
                  s1, signif(t1, 5),
                  p1, stars),
                  ncol = 6)
  colnames(obj) <- c(" ", "Estimate", "Std.Error", "t-value", "Pr(>|t|)", " ")
  obj <- as.data.frame(obj)

  ## models
  w <- object$modelSummary$w
  nMod <- length(w)
  ylags <- object$ylags
  if(is.null(ylags)){
    ylags <- matrix(0, ncol = nMod)
  }else if(is.vector(ylags)){
    ylags <- matrix(ylags, ncol = 1)
  }
  mulags <- object$mulags
  if(is.null(mulags)){
    mulags <- matrix(0, ncol = nMod)
  }else if(is.vector(mulags)){
    mulags <- matrix(mulags, ncol = 1)
  }
  r <- object$modelSummary$r

  for(i in seq_len(nrow(object$modelSummary))){
    if(nrow(object$modelSummary) > 1){
      cat(paste0("Model ", i, ":\n"))
    }else{
      cat("Model:\n")
    }
    cat(paste("intercept:            ", w[i], "\n"))
    cat(paste("observation lags:     ", paste(ylags[,i], collapse = ", "), "\n"))
    cat(paste("mean lags:            ", paste(mulags[,i], collapse = ", "), "\n"))
    cat(paste("exogenous covariates: ", r[i], "\n"))
    cat("\n")
  }

  ## print coefficient matrix
  cat("Coefficients: \n")
  print(obj, print.gap = 2L, digits = 5, quote = FALSE, row.names = FALSE)
  cat("----- \nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  cat(paste("\nLagged parameters tested on t-distribution with", n1, "degrees of freedom"))
  cat("\nWarning: t-distribution may not be a good approximation due to constraints")
  if(exists("n2"))
    cat(paste("\nRho tested on t-distribution with", n2, "degrees of freedom"))
  cat("\n----- \n \n")

  ## Number of observations, etc.
  cat(paste("Number of: \n  Observations : ", signif(object$nobs, 5),
            "\n  Time Series  : ", signif(object$nts, 5),
            "\n  Parameters   : ", signif(object$npar, 5)))
  cat("\n")
  cat(paste("\nLog likelihood : ", signif(object$logLik, 5),
            "\nAIC            : ", signif(2*(object$npar - object$logLik), 5),
            "\nBIC            : ", signif(object$npar*log(object$nobs) - 2*object$logLik, 5)))
  cat("\n")
  invisible(object)
}

#' @rdname PoARX_methods
#'
#' @export
vcov.PoARX <- function(object, ...){
  if(is.null(object$varCovar)){
    message("Standard errors not available")
    NULL
  }else{
    object$varCovar
  }
}

#' @title Predict method for PoARX objects
#'
#' @description The \code{predict} function can be used to predict the next h steps
#'    given new data or return fitted value. If \code{newdata} is \code{NULL}, then
#'    fitted values will be returned if \code{"response"} is chosen, otherwise
#'    the probability of each observed data point. Missing values can be represented by
#'    \code{NA} and are dealt with by \code{PoARXPredict}. Note that as h increases,
#'    the variances of the predictions will increase. Confidence intervals can be
#'    returned if desired, using simulation.
#'
#' @param object A \code{"PoARX"} object from \code{PoARXFit}.
#' @param newdata An optional data.frame for the new covariates and observations.
#'    The names of the data.frame should match the original variable names. If some new
#'    observations have not been obtained, put an NA value where the missing value is
#'    and the observation will be estimated by \code{PoARXPredict}.
#' @param type A character argument indicating the type of response. If \code{"response"},
#'    return the mean parameter for each time series. If \code{"probability"}, give the
#'    probability of the response.
# #' @param interval A logical argument for the calculation of confidence intervals.
# #' @param cl A numeric value representing the desired confidence level for predictions of
# #'    response values.
# #' @param B A number of desired simulations in the calculation of confidence intervals in
# #'    any \code{'response'} values.
# #' @param trace A logical indicator for the printing of update messages.
#' @param ... Further arguments.
#'
#' @return If \code{type = "response"} a data.frame containing
#'    \eqn{\lambda} values for each time series, along with optional confidence intervals.
#'
#'    If \code{type = "probability"} and the respective observations are provided, then
#'    a matrix of probabilities is returned. If no observations are given, there
#'    are too many possible values, so the \eqn{\lambda}{lambda} values are returned instead,
#'    as if \code{type = "response"}.
#'
#'
#' @examples
#' #### Examples go here!
#'
#' @export
predict.PoARX <- function(object, newdata = NULL, type = c("response", "probability"),
                          # interval = FALSE, cl = 0.95, B = 1000, trace = TRUE,
                          ...){

  #### I've suspended the interval bit due to comments in JTSA doc

  type <- match.arg(type)
  # if(interval & type == "probability")
  #   stop("intervals are only available for 'response' predictions")

  if(missing(newdata) || is.null(newdata)){
    #### If no new data
    switch(type,
           response = {
             #### Return fitted values
           #  if(interval){
          #     olddata <- cbind(object$Y, object$X)
               # res <- .predictLambda(object, olddata, interval, cl, B, trace)
          #     res <- .predictLambda(object, olddata)
          #   }else{
               res <- object$fitted.values
           #  }
           },
           probability = {
             res <- object$fitted.probs
           },
           stop("'type' must be one of 'response' and 'probability'")
           )
    }else{
        switch(type,
               response = {
                 # newlam <- .predictLambda(object, newdata, interval, cl, B, trace)
                 newlam <- .predictLambda(object, newdata)
                 res <- newlam
               },
               probability = {
                 newlam <- .predictLambda(object, newdata)#, trace = FALSE)
                 newlam <- as.matrix(newlam)

                 Ynew <- NULL
                 if(all(names(object$Y) %in% names(newdata))){
                   Ynew <- model.part(object$formula, data = newdata,
                                      lhs = seq_len(object$nts))
                   Ynew <- as.matrix(Ynew)
                 }
                 if(is.null(Ynew)){
                   message("Cannot calculate probabilities - too many possibilities.\nFitted values given.")
                   res <- newlam
                 }else{
                   ind <- !is.na(Ynew)
                   ind <- apply(ind, 1, function(x) all(x))
                   prob <- matrix(NA, nrow(Ynew), 1)
                   if(object$rhoInd){
                     prob[ind,] <- dFranksCopulaPoisson(Ynew[ind , , drop = FALSE],
                                                        newlam[ind , , drop = FALSE],
                                                        object$rhoValue)
                   }else{
                     obsProb <- dpois(Ynew[ind , , drop = FALSE],
                                      newlam[ind , , drop = FALSE])
                     obsProb <- apply(obsProb, 1, prod)
                     prob[ind,] <- obsProb
                   }
                   cond <- any(!ind)
                   if(cond)
                     warning("probabilities after NA are calculated based on observed value equalling intensity!")
                   res <- prob
                 }
               },
               stop("'type' must be one of 'response' or 'probability'")
               )
    }
  res
}
