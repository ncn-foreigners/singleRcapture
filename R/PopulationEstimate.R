#' Total Population Point and Interval estimate
#'
#' Creates a Horvitz-Thompson/Chao/Zelterman
#' point and interval estimate for total
#' and missing population from zero truncated poisson model.
#'
#' @param y Observed values
#' @param X A matrix of covariates
#' @param Grad A gradient of a model with respect to regression parameters
#' @param Parameter An estimated Parameter for the model
#' @param beta Fitted regression values
#' @param family Model Family
#' @param weights If model is weighted weights of particular observations
#' @param Hess Hessian of a model
#' @param dispersion Estimated dispersion parameter
#' for truncated Negative binomial distributions
#' @param Method A method of constructing confidence interval either analytic
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
#' @export
#'
PopulationEstimate <- function(y,
                               X,
                               Grad,
                               Parameter,
                               beta,
                               weights = 1,
                               Hess,
                               family,
                               trcount,
                               dispersion,
                               Method = "analytic") {
  if (Method == "analytic") {
    N <- family$Point.est(disp = dispersion,
                          pw = weights,
                          Lambda = Parameter) + trcount

    Variation <- family$Pop.var(beta = beta, pw = weights,
                                Lambda = Parameter,
                                disp = dispersion,
                                Hess = Hess, X = X)

    Confidence_Interval <- c(Lower_Bound = max(N - 1.96 * sqrt(Variation),
                                              (length(y) + trcount)),
                             Upper_Bound = N + 1.96 * sqrt(Variation))
  } else if (grepl("bootstrap", Method, fixed = TRUE)) {
    if (family$family %in% c("chao", "zelterman")) {
      N <- family$Point.est(disp = dispersion,
                            pw = weights,
                            Lambda = Parameter) + trcount
      f0 <- N - sum(weights) - trcount
    } else {
      N <- family$Point.est(disp = dispersion,
                            pw = weights,
                            Lambda = Parameter) + trcount
      f0 <- N - length(y) - trcount
    }
    if (trcount == 0 || is.null(trcount)) {
      Prob <- c(f0, as.numeric(table(y))) / N
      Prob <- cumsum(Prob)
      names(Prob) <- c(0, as.numeric(names(table(y))))
    } else {
      if (family$family %in% c("chao",
                               "zelterman")) {
        Prob <- c(f0, weights, trcount) / N
        Prob <- cumsum(Prob)
        names(Prob) <- c(0, as.numeric(names(table(y))), 3)
      } else {
        Prob <- c(f0, trcount, as.numeric(table(y))) / N
        Prob <- cumsum(Prob)
        names(Prob) <- c(0, 1, as.numeric(names(table(y))))
      }
    }
    StrappedStatistic <- NULL

    Bootnumber <- 1000

    for (z in 1:Bootnumber) {
      U <- stats::runif(N)
      Strap <- NULL

      for (k in names(Prob)) {
        Strap <- c(Strap, length(which(U <= as.numeric(Prob[k]))))
      }
      names(Strap) <- names(Prob)

      m <- length(Strap)

      while (m != 1) {
        Strap[m] <- Strap[m] - Strap[m - 1]
        m <- m - 1
      }
      full <- NULL

      for (k in names(Strap)) {
        full <- c(full, rep(as.numeric(k), Strap[k]))
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
        Theta <- stats::optimize(f = ll,
                                 interval = c(-start,2 * start))$minimum
        Theta <- family$linkinv(matrix(c(1, 1), ncol = 1) %*% Theta)
      } else{
        Theta <- IRLS(Dependent = as.numeric(full),
                      Covariates = matrix(rep(1, length(full)), ncol = 1),
                      family = family,
                      start = start,
                      disp = dispersion,
                      disp.given = TRUE)$Coefficients
        Theta <- family$linkinv(matrix(rep(1, length(full)), ncol = 1) %*% Theta)
      }

      if (family$family %in% c("chao", "zelterman")) {
        StrappedStatistic <- c(StrappedStatistic,
                               family$Point.est(disp = dispersion,
                                                pw = as.numeric(table(full)),
                                                Lambda = Theta) + trcountboot)
      } else {
        StrappedStatistic <- c(StrappedStatistic,
                               family$Point.est(disp = dispersion,
                                                pw = 1,
                                                Lambda = Theta) + trcountboot)
      }
    }
    N <- mean(StrappedStatistic)
    if (Method == "bootstrapSD") {
      Variation <- stats::var(StrappedStatistic)
      Confidence_Interval <- c(Lower_Bound = max(N - 1.96 * sqrt(Variation),
                                                 (length(y) + trcount)),
                               Upper_Bound = N + 1.96 * sqrt(Variation))
    } else if (Method == "bootstrapPerc") {
      Variation <- stats::var(StrappedStatistic)
      Confidence_Interval <- stats::quantile(StrappedStatistic, c(0.025, 0.975))
      names(Confidence_Interval) <- c("Lower_Bound", "Upper_Bound")
    }
  }

  list(Point_estimate = N,
       Variance = Variation,
       Confidence_Interval = Confidence_Interval)
}
