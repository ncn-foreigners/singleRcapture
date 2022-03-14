#' Total Population Point and Interval estimate
#'
#' Creates a Horvitz-Thompson/Chao/Zelterman
#' point and interval estimate for total
#' and missing population from zero truncated poisson model.
#'
#' @param y Observed values
#' @param X A matrix of covariates
#' @param grad A gradient of a model with respect to regression parameters
#' @param parameter An estimated Parameter for the model
#' @param beta Fitted regression values
#' @param family Model Family
#' @param weights If model is weighted weights of particular observations
#' @param hessian Hessian of a model
#' @param dispersion Estimated dispersion parameter
#' for truncated Negative binomial distributions
#' @param method A method of constructing confidence interval either analytic
#' to use formula for analytic CI or bootstrap where bootstraped confidence
#' interval may either be based on 2.5%-97.5% percientiles ("bootstrapPerc")
#' or by estimating SD ("bootstrapSD")
#' @param trcount Optional parameter for Zero-one truncated models, if population estimate
#' is for one inflated model then it specifies one counts and includes them in
#' final population estimate both point and interval, and for zeltermann/chao
#' estimator where it specifies counts of not used in estimate
#'
#' @return Returns a list of size 3 with:
#' Point estimate, Interval estimate, and Variance
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom stats quantile
#' @export
#'
populationEstimate <- function(y,
                               X,
                               grad,
                               parameter,
                               beta,
                               weights = 1,
                               hessian,
                               family,
                               trcount,
                               dispersion,
                               method = "analytic") {
  if (method == "analytic") {
    N <- family$pointEst(disp = dispersion,
                          pw = weights,
                          lambda = parameter) + trcount

    variation <- family$popVar(beta = beta, pw = weights,
                                lambda = parameter,
                                disp = dispersion,
                                hess = hessian, X = X)

    confidenceInterval <- c(lowerBound = max(N - 1.96 * sqrt(variation),
                                            (length(y) + trcount)),
                             upperBound = N + 1.96 * sqrt(variation))
  } else if (grepl("bootstrap", method, fixed = TRUE)) {

    if (family$family %in% c("chao", "zelterman")) {
      N <- family$pointEst(disp = dispersion,
                           pw = weights,
                           lambda = parameter) + trcount
      f0 <- N - sum(weights) - trcount
    } else {
      N <- family$pointEst(disp = dispersion,
                           pw = weights,
                           lambda = parameter) + trcount
      f0 <- N - length(y) - trcount
    }
    if (trcount == 0 || is.null(trcount)) {
      prob <- c(f0, as.numeric(table(y))) / N
      prob <- cumsum(prob)
      names(prob) <- c(0, as.numeric(names(table(y))))
    } else {
      if (family$family %in% c("chao",
                               "zelterman")) {
        prob <- c(f0, weights, trcount) / N
        prob <- cumsum(prob)
        names(prob) <- c(0, as.numeric(names(table(y))), 3)
      } else {
        prob <- c(f0, trcount, as.numeric(table(y))) / N
        prob <- cumsum(prob)
        names(prob) <- c(0, 1, as.numeric(names(table(y))))
      }
    }
    strappedStatistic <- NULL

    bootnumber <- 10000

    for (z in 1:bootnumber) {
      U <- stats::runif(N)
      strap <- NULL

      for (k in names(prob)) {
        strap <- c(strap, length(which(U <= as.numeric(prob[k]))))
      }
      names(strap) <- names(prob)

      m <- length(strap)

      while (m != 1) {
        strap[m] <- strap[m] - strap[m - 1]
        m <- m - 1
      }
      full <- NULL

      for (k in names(strap)) {
        full <- c(full, rep(as.numeric(k), strap[k]))
      }

      full <- full[full != 0]
      trcountboot <- 0
      if (family$family %in% c("zotpoisson", "zotnegbin")) {
        trcountboot <- length(full[full == 1])
        full <- full[full != 1]
      } else if (family$family %in% c("chao", "zelterman")) {
        trcountboot <- length(full[full == 3])
        full <- full[full < 3]
      }

      start <- beta
      if (!is.null(dispersion)) {
        start <- beta[-1]
      }

      start <- mean(start)

      if (family$family %in% c("chao", "zelterman")) {
        ll <- family$make_minusloglike(y = c(1, 2),
                                       X = matrix(c(1, 1), ncol = 1),
                                       weight = as.numeric(table(full)))
        gr <- family$make_gradient(y = c(1, 2),
                                   X = matrix(c(1, 1), ncol = 1),
                                   weight = as.numeric(table(full)))
        Theta <- stats::optim(par = start,
                              lower = -start,
                              upper = 2 * start,
                              fn = ll,
                              gr = function (x) -gr(x),
                              method = "Brent",
                              control = list(reltol = .Machine$double.eps))$par
        Theta <- family$linkinv(matrix(c(1, 1), ncol = 1) %*% Theta)
      } else{
        theta <- IRLS(dependent = as.numeric(full),
                      covariates = matrix(rep(1, length(full)), ncol = 1),
                      family = family,
                      start = start,
                      disp = dispersion,
                      disp.given = TRUE)$coefficients
        theta <- family$linkinv(matrix(rep(1, length(full)), ncol = 1) %*% theta)
      }

      if (family$family %in% c("chao", "zelterman")) {
        strappedStatistic <- c(strappedStatistic,
                               family$pointEst(disp = dispersion,
                                               pw = as.numeric(table(full)),
                                               lambda = theta) + trcountboot)
      } else {
        strappedStatistic <- c(strappedStatistic,
                               family$pointEst(disp = dispersion,
                                               pw = 1,
                                               lambda = theta) + trcountboot)
      }
    }

    if (method == "bootstrapPerc") {
      variation <- stats::var(strappedStatistic)
      confidenceInterval <- stats::quantile(strappedStatistic,
                                            c(0.025, 0.975))
      names(confidenceInterval) <- c("lowerBound", "upperBound")

    } else {
      variation <- stats::var(strappedStatistic)
      confidenceInterval <- c(lowerBound = max(N - 1.96 * sqrt(variation),
                                               (length(y) + trcount)),
                              upperBound = N + 1.96 * sqrt(variation))
    }
  }

  list(pointEstimate = N,
       variance = variation,
       confidenceInterval = confidenceInterval)
}
