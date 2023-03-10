#' @title Find the constraints for PoARX model(s) of the same structure
#'
#' @description
#' There are various constraints placed on the parameters \eqn{\omega},
#'    {\eqn{\alpha_i}}, {\eqn{\beta_i}}, and \eqn{\eta} in a PoARX model.
#'    This function creates matrices which can be used to implement such
#'    constraints in \code{constrOptim}.
#'
#' @param M the number of sets of PoARX coefficients required.
#' @param r the number of exogenous covariates.
#' @param rho a logical indicator for whether a column is needed for a dependence
#'    parameter.
#' @param intercept A logical indicator for whether an intercept should be
#'    included in the model.
#' @param ylags A vector of time lags in the filter equation relating to the
#'    observations.
#' @param mulags A vector of time lags in the filter equation relating to the
#'    intensities.
#' @param ... currently not used.
#'
#' @return A list containing the following:
#'     \item{ui}{A constraint matrix.}
#'     \item{ci}{A constraint vector.}
#'
#' @details
#' The conditions for the parameters in a PoARX model are as follows. Firstly,
#'    all parameters must be non-negative to ensure that the Poisson parameter
#'    \eqn{\lambda} is non-negative. Secondly, the autoregressive parameters
#'    {\eqn{\alpha_i}} and {\eqn{\beta_i}} must satisfy
#'       \deqn{ \sum (\alpha_i + \beta_i) < 1 }
#'    in order to abide by the stationarity and ergodicity conditions.
#'
#' @export
#'
PoARXConstraints <- function(M, r, rho = FALSE, intercept = TRUE,
                             ylags = NULL, mulags = NULL, ...){
  w <- as.numeric(intercept)
  P <- length(ylags)
  Q <- length(mulags)
  N <- w + P + Q + r

  ui <- diag(1, nrow = N*M, ncol = N*M)
  ci <- rep(0, N*M)
  if(P + Q > 0){
    ui <- rbind(ui, matrix(0, nrow = M, ncol = N*M))
    for(j in seq_len(M)){
      ui[N*M + j, (N*(j-1) + w + min(1, P + Q)):(N*(j-1) + w + P + Q)] <- -1
    }
    ci <- c(ci, rep(-1, M))
  }
  if(rho){
    ui <- cbind(ui, 0)
    if(M>2){
      ui <- rbind(ui, 0)
      ui[nrow(ui), ncol(ui)] <- 1
      ci <- c(ci, 0)
    }
  }

  list(ui = ui, ci = ci)
}
