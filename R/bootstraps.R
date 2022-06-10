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
#' @param trace TODO
#' @param control.bootstrap.method TODO
#' @param method TODO
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
                      lambda,
                      trace,
                      control.bootstrap.method = NULL,
                      method) {
  strappedStatistic <- NULL
  n <- length(y)
  famName <- family$family
  weights <- weights[[1]]
  
  if (length(weights) == 1) {
    weights <- rep(1, length(lambda))
  }
  
  for (k in 1:numboot) {
    strap <- sample.int(replace = TRUE, n = n)
    ystrap <- as.numeric(y[strap])
    weightsstrap <- as.numeric(weights[strap])
    Xstrap <- as.matrix(X[strap, ])
    
    weightsstraptemp <- weightsstrap
    Xstraptemp <- Xstrap
    
    if (famName == "zelterman") {
      df <- data.frame(ystrap,
                       weightsstrap,
                       Xstrap)

      df <- df[df[1] == 1 | df[1] == 2, ]

      ystrap <- as.numeric(unlist(df[1]))
      weightsstraptemp <- as.numeric(unlist(df[2]))
      Xstraptemp <- as.matrix(df[-c(1, 2)])
    } else if (grepl(x = famName, pattern = "^zot.*")) {
      df <- data.frame(ystrap,
                       weightsstrap,
                       Xstrap)
      
      trcount <- nrow(df[df["ystrap"] == 1, ])
      df <- df[df["ystrap"] > 1, ]
      ystrap <- as.numeric(unlist(df[1]))
      weightsstraptemp <- as.numeric(unlist(df[2]))
      Xstraptemp <- as.matrix(df[-c(1, 2)])
    } else if (famName == "chao") {
      df <- data.frame(ystrap,
                       weightsstrap,
                       Xstrap)
      
      df <- df[df[1] == 1 | df[1] == 2, ]
      
      ystrap <- as.numeric(df[, 1])
      weightsstraptemp <- as.numeric(unlist(df[, 2]))
      weightsstrap <- weightsstraptemp
      Xstraptemp <- as.matrix(df[, -c(1, 2)])
      Xstrap <- Xstraptemp
    }
    
    theta <- NULL
    try(
      theta <- estimate_popsize.fit(
        y = ystrap,
        X = Xstraptemp,
        family = family,
        start = beta,
        dispersion = dispersion,
        method = method,
        prior.weights = weightsstraptemp,
        control = control.bootstrap.method
      )$beta,
      silent = TRUE
    )
    
    if (is.null(theta)) {
      k <- k - 1
    } else {
      theta <- family$linkinv(Xstrap %*% theta)
      if (isTRUE(trace)) {print(theta)}
      
      strappedStatistic <- c(strappedStatistic,
                             family$pointEst(disp = dispersion,
                                             pw = weightsstrap,
                                             lambda = theta) + trcount)
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
#' @param trace TODO
#' @param control.bootstrap.method TODO
#' @param method TODO
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
                    lambda,
                    trace,
                    control.bootstrap.method = NULL,
                    method) {
  strappedStatistic <- NULL
  n <- length(y)
  famName <- family$family
  if (length(weights[[1]]) == 1) {
    weights[[1]] <- rep(1, length(lambda))
  }
  
  if (length(weights[[2]]) == 1) {
    weights[[2]] <- rep(1, n)
  }
  
  if (family$family %in% c("chao", "zelterman")) {
    cat("Probability model will be taken as given by poisson distribution",
    "since zelterman and chao models are based on mixture of poisson",
    "distribution semi-parametric bootstrap may be a better choice.", 
    sep = "\n")
  } else if (family$family %in% c("ztnegbin",
                                  "zotnegbin")) {
    cat("Due to many possible computational problems with truncated binomial models",
    "bootstrap samples will be drawn until specified number of them can be fitted", 
    "this may significantly contribute to increase in runtime.", 
    sep = "\n")
  }
  dataFunc <- switch(famName,
  "ztpoisson"  = function(lambda, n, disp) {stats::rpois(n = n, lambda = lambda)},
  "chao"       = function(lambda, n, disp) {stats::rpois(n = n, lambda = lambda)},
  "zelterman"  = function(lambda, n, disp) {stats::rpois(n = n, lambda = lambda)},
  "zotpoisson" = function(lambda, n, disp) {stats::rpois(n = n, lambda = lambda)},
  "ztnegbin"   = function(lambda, n, disp) {stats::rnbinom(n = n, mu = lambda, size = exp(-disp))},
  "zotnegbin"  = function(lambda, n, disp) {stats::rnbinom(n = n, mu = lambda, size = exp(-disp))},
  "ztgeom"     = function(lambda, n, disp) {stats::rgeom(n = n, prob = (1 / (1 + lambda)))},
  "zotgeom"    = function(lambda, n, disp) {stats::rgeom(n = n, prob = (1 / (1 + lambda)))})

  if (isFALSE(grepl(x = famName, pattern = "^(zot|cha).*"))) {
  contr <- family$pointEst(disp = dispersion, 
                           pw = weights[[2]],
                           lambda = lambda,
                           contr = TRUE)
  } else {
    contr <- vector(mode = "numeric", length = length(y))
    cond <- switch(
      substr(famName, 1, 3),
      "zot" = (y != 1),
      "cha" = (y < 3)
    )
    contr[cond] <- family$pointEst(disp = dispersion, 
                                   pw = weights[[1]], 
                                   lambda = lambda,
                                   contr = TRUE)
    contr[!cond] <- 1
  }
  N <- round(sum(contr))
  prob <- contr - floor(contr)
  contr <- floor(contr) + stats::rbinom(n = n, size = 1, prob = prob)
  
  prob <- contr / N
  weights <- weights[[2]]
  getlambda <- family$linkinv

  for (k in 1:numboot) {
    strap <- sample.int(replace = TRUE, n = n, size = N, prob = prob)
    weightsstrap <- as.numeric(weights[strap])
    Xstrap <- as.matrix(X[strap, ])

    #print(summary(family$linkinv(Xstrap %*% beta)))
    ystrap <- dataFunc(n = N,
                       lambda = getlambda(Xstrap %*% beta),
                       disp = dispersion)
    weightsstraptemp <- weightsstrap
    Xstraptemp <- Xstrap

    df <- data.frame(ystrap,
                     weightsstraptemp,
                     Xstraptemp)
    trcount <- 0
    if (grepl(x = famName, pattern = "^zot.*")) trcount <- nrow(df[df["ystrap"] == 1, ])
    ifelse(grepl(x = famName, pattern = "^zot.*"),
           df <- df[df["ystrap"] > 1, ],
           df <- df[df["ystrap"] > 0, ])

    if (isTRUE(trace)) cat("Iteration number:", k, "untruncated sample size:", nrow(df), "\n", sep = " ")
    ystrap <- as.numeric(unlist(df[, 1]))
    weightsstrap <- weightsstraptemp <- as.numeric(unlist(df[, 2]))
    Xstrap <- Xstraptemp <- as.matrix(df[, -c(1, 2)])
    
    if (famName == "zelterman") {
      df <- data.frame(ystrap,
                       weightsstrap,
                       Xstrap)
      
      df <- df[df[1] == 1 | df[1] == 2, ]
      
      ystrap <- as.numeric(df[, 1])
      weightsstraptemp <- as.numeric(unlist(df[, 2]))
      Xstraptemp <- as.matrix(df[, -c(1, 2)])
    } else if (famName == "chao") {
      df <- data.frame(ystrap,
                       weightsstrap,
                       Xstrap)
      
      trcount <- nrow(df[df[1] > 2, ])

      df <- df[df[1] == 1 | df[1] == 2, ]
      
      ystrap <- as.numeric(df[, 1])
      weightsstraptemp <- as.numeric(unlist(df[, 2]))
      weightsstrap <- weightsstraptemp
      Xstraptemp <- as.matrix(df[, -c(1, 2)])
      Xstrap <- Xstraptemp
    }
    
    theta <- NULL
    try(
      if (method == "mle") {
        methodopt <- control.bootstrap.method$mleMethod
        log_like <- family$make_minusloglike(y = ystrap, X = Xstraptemp, weight = weightsstraptemp)
        grad <- family$make_gradient(y = ystrap, X = Xstraptemp, weight = weightsstraptemp)
        ctrl <- control.bootstrap.method$optimPass
        if (is.null(ctrl)) {
          list(maxit = control.bootstrap.method$maxiter,
               factr = control.bootstrap.method$epsilon)
        }
        theta <- stats::optim(
          par = if (grepl(x = family$family, pattern = "negbin")) c(dispersion, beta) else beta,
          fn = log_like,
          gr = function(x) -grad(x),
          method = methodopt,
          control = ctrl
        )$par
        
        if (grepl(x = family$family, pattern = "negbin")) {theta <- theta[-1]}
      } else if (method == "robust") {
        theta <- IRLS(
          dependent = ystrap,
          family = family,
          covariates = Xstraptemp,
          start = beta,
          disp.given = TRUE,
          disp = dispersion,
          eps = control.bootstrap.method$epsilon,
          maxiter = control.bootstrap.method$maxiter,
          silent = TRUE,
          trace = FALSE
        )$coefficients
      },
      silent = TRUE
    )

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
#' @param trace TODO
#' @param control.bootstrap.method TODO
#' @param method TODO
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
                       lambda,
                       trace,
                       control.bootstrap.method = NULL,
                       method) {
  strappedStatistic <- NULL
  n <- length(y)
  weights <- weights
  stop("Not yet implemented")
  if (length(weights) == 1) {
    weights <- rep(1, n)
  }
  
  for (k in 1:numboot) {
    
  }
  
  strappedStatistic
}