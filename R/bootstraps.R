#' Bootstraps for population size variation estimation
#'
#' @description {
#'  TODO
#' }
#' @param family TODO
#' @param y TODO
#' @param X TODO
#' @param dispersion TODO
#' @param beta TODO
#' @param weights TODO
#' @param trcount TODO
#' @param numboot TODO
#'
#' @return TODO
#' @export
noparBoot <- function(family,
                      y,
                      X,
                      dispersion,
                      beta,
                      weights,
                      trcount,
                      numboot) {
  strappedStatistic <- NULL
  n <- length(y)
  famName <- family$family
  
  if (length(weights) == 1) {
    weights <- rep(1, n)
  }
  
  for (k in 1:numboot) {
    strap <- sample.int(replace = TRUE, n = n)
    ystrap <- as.numeric(y[strap])
    weightsstrap <- as.numeric(weights[strap])
    Xstrap <- as.matrix(X[strap, ])
    
    ystraptemp <- ystrap
    weightsstraptemp <- weightsstrap
    Xstraptemp <- Xstrap
    
    if (famName == "zelterman") {
      df <- data.frame(ystrap,
                       weightsstrap,
                       Xstrap)

      df <- df[df[1] == 1 | df[1] == 2, ]

      ystraptemp <- as.numeric(unlist(df[1]))
      weightsstraptemp <- as.numeric(unlist(df[2]))
      Xstraptemp <- as.matrix(df[-c(1, 2)])
    }
    
    theta <- IRLS(dependent = ystraptemp,
                  covariates = Xstraptemp,
                  family = family,
                  start = beta,
                  disp = dispersion,
                  disp.given = TRUE,
                  eps = 1e-6,
                  weights = weightsstraptemp,
                  maxiter = 50,
                  silent = TRUE)$coefficients
    
    theta <- family$linkinv(Xstrap %*% theta)
    
    strappedStatistic <- c(strappedStatistic,
                           family$pointEst(disp = dispersion,
                                           pw = weightsstrap,
                                           lambda = theta) + trcount)
  }
  
  strappedStatistic
}
#' Bootstraps for population size variation estimation
#'
#' @description {
#'  TODO
#' }
#' @param family TODO
#' @param y TODO
#' @param X TODO
#' @param dispersion TODO
#' @param beta TODO
#' @param weights TODO
#' @param trcount TODO
#' @param numboot TODO
#'
#' @return TODO
#' @export
parBoot <- function(family,
                    y,
                    X,
                    dispersion,
                    beta,
                    weights,
                    trcount,
                    numboot) {
  strappedStatistic <- NULL
  n <- length(y)
  if (length(weights) == 1) {
    weights <- rep(1, n)
  }
  
  for (k in 1:numboot) {
    
  }
  
  strappedStatistic
}