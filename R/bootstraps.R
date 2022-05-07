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
#' @param lambda TODO
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
                      numboot,
                      lambda) {
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
#' @param lambda TODO
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
                    numboot,
                    lambda) {
  strappedStatistic <- NULL
  n <- length(y)
  famName <- family$family
  
  if (length(weights) == 1) {
    weights <- rep(1, n)
  }
  
  if (family$family %in% c("chao", "zelterman")) {
    warning("Probability model will be taken as given by poisson distribution,\nsince zelterman and chao models are based on mixture of poisson\ndistribution semi-parametric bootstrap may be a better choice.")
  } else if (family$family %in% c("ztnegbin",
                                  "zotnegbin")) {
    cat("Due to many possible computational problems with truncated binomial models",
    "bootstrap samples will be drawn until specified number of them can be fitted", 
    "this may significantly contribute to increase in runtime.", sep = "\n")
  }
  
  contr <- family$pointEst(disp = dispersion, 
                           pw = weights,
                           lambda = lambda,
                           contr = TRUE)
  N <- round(sum(contr))
  
  prob <- contr - floor(contr)
  contr <- floor(contr) + stats::rbinom(n = n, size = 1,
                                        prob = prob)
  dataFunc <- ifelse(famName %in% c("ztpoisson", 
                                    "chao",
                                    "zelterman",
                                    "zotpoisson"),
                     function(lambda, n, disp) {stats::rpois(n = n, 
                                                             lambda = lambda)},
                     ifelse(famName %in% c("ztnegbin", 
                                           "zotnegbin"),
                            function(lambda, n, disp) {stats::rnbinom(n = n, 
                                                                      mu = lambda,
                                                                      size = exp(-disp))},
                            ifelse(famName %in% c("ztgeom",
                                                  "zotgeom"),
                                   function(lambda, n, disp) {stats::rgeom(n = n,
                                                                           prob = (1 / (1 + lambda)))},
                                   "")))
    
  prob <- contr / N
  
  for (k in 1:numboot) {
    strap <- sample.int(replace = TRUE, n = n, size = N, prob = prob)
    weightsstrap <- as.numeric(weights[strap])
    Xstrap <- as.matrix(X[strap, ])
    
    ystraptemp <- dataFunc(n = N,
                           lambda = exp(Xstrap %*% beta),
                           disp = dispersion)
    weightsstraptemp <- weightsstrap
    Xstraptemp <- Xstrap
    
    df <- data.frame(ystraptemp,
                     weightsstraptemp,
                     Xstraptemp)
    ifelse(famName %in% c("zotpoisson",
                          "zotgeom",
                          "zotnegbin"),
           df <- df[df["ystraptemp"] > 1, ],
           df <- df[df["ystraptemp"] > 0, ])
    
    ystraptemp <- as.numeric(df[, 1])
    ystrap <- ystraptemp
    weightsstraptemp <- as.numeric(unlist(df[, 2]))
    weightsstrap <- weightsstraptemp
    Xstraptemp <- as.matrix(df[, -c(1, 2)])
    Xstrap <- Xstraptemp
    
    if (famName == "zelterman") {
      df <- data.frame(ystrap,
                       weightsstrap,
                       Xstrap)
      
      df <- df[df[1] == 1 | df[1] == 2, ]
      
      ystraptemp <- as.numeric(df[, 1])
      weightsstraptemp <- as.numeric(unlist(df[, 2]))
      Xstraptemp <- as.matrix(df[, -c(1, 2)])
    } else if (famName == "chao") {
      df <- data.frame(ystrap,
                       weightsstrap,
                       Xstrap)
      
      df <- df[df[1] == 1 | df[1] == 2, ]
      
      ystraptemp <- as.numeric(df[, 1])
      ystrap <- ystraptemp
      weightsstraptemp <- as.numeric(unlist(df[, 2]))
      weightsstrap <- weightsstraptemp
      Xstraptemp <- as.matrix(df[, -c(1, 2)])
      Xstrap <- Xstraptemp
    }

    ll <- family$make_minusloglike(y = ystraptemp,
                                   X = Xstraptemp,
                                   weight = weightsstraptemp)
    gr <- family$make_gradient(y = ystraptemp,
                               X = Xstraptemp,
                               weight = weightsstraptemp)
    
    if (family$family %in% c("ztnegbin",
                             "zotnegbin")) {
      theta <- NULL
      try(theta <- stats::optim(par = c(dispersion, beta),
                                fn = ll,
                                gr = function(x) -gr(x),
                                control = list(reltol = 1e-5,
                                               warn.1d.NelderMead = FALSE,
                                               maxit = 50))$par[-1],
          silent = TRUE)
    } else if (family$family %in% c("ztgeom",
                                    "zotgeom")) {
      theta <- stats::optim(par = beta,
                            fn = ll,
                            gr = function(x) -gr(x),
                            control = list(reltol = 1e-5,
                                           warn.1d.NelderMead = FALSE,
                                           maxit = 50))$par
    } else {
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
    }
    
    if(!is.null(theta)) {
      theta <- family$linkinv(Xstrap %*% theta)
      
      strappedStatistic <- c(strappedStatistic,
                             family$pointEst(disp = dispersion,
                                             pw = weightsstrap,
                                             lambda = theta) + trcount)
    } else {
      k <- k - 1
    }
    
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
#' @param lambda TODO
#'
#' @return TODO
#' @export
semparBoot <- function(family,
                       y,
                       X,
                       dispersion,
                       beta,
                       weights,
                       trcount,
                       numboot,
                       lambda) {
  strappedStatistic <- NULL
  n <- length(y)
  if (length(weights) == 1) {
    weights <- rep(1, n)
  }
  
  for (k in 1:numboot) {
    
  }
  
  strappedStatistic
}