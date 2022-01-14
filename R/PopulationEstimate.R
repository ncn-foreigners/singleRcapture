#' Total Population Point and Interval estimate
#'
#' Creates a Horwitz-Thompson point and interval estimate for total and missing population from zero truncated poisson model.
#'
#' @param y Observed values
#' @param X A covariate matrix
#' @param Grad A gradient of a model with respect to regression parameters
#' @param Parameter An estimated Parameter for the model
#' @param beta Fitted regression values
#' @param family Model Family
#' @param weights If model is weighted weights of particular observations
#' @param Hess Hessian of a model
#' @param Overdispersion For truncated Negative binomial distributions
#' @param Method A method of constructing confidence interval either analytic or bootstrap where bootstraped confidence interval may
#' either be based on 2.5%-97.5% percientiles ("bootstrapPerc") or by estimating SD ("bootstrapSD") to use formula for analytic CI
#'
#' @return Returns a list of size 3 with: \cr
#' (1) Point estimate \cr
#' (2) Variance \cr
#' (3) Confidence interval constructed using Point estimate plus minus 1.96
#' (approximate of 97.5% percentile of standardaised normal distribution) times standard deviation
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
                               Overdispersion = NULL,
                               Method = "analytic") {
  if (family$family == "ZTP"){
    Lambda <- Parameter
    N <- sum(1 / (1 - exp(-Lambda)))

    if (Method == "analytic") {
      X <- as.data.frame(X)
      Inform <- -Hess(beta)

      f1 <- colSums(-X * (exp(log(Lambda)-Lambda)/((1 - exp(-Lambda)) ** 2)))
      f1 <- t(f1) %*% solve(as.matrix(Inform)) %*% f1

      f2 <- sum(exp(-Lambda)/((1-exp(-Lambda)) ** 2))

      Variation <- f1 + f2

      Confidence_Interval <- c(Lower_Bound = max(N - 1.96 * sqrt(Variation), length(y)),
                               Upper_Bound = N + 1.96 * sqrt(Variation))
    } else if (grepl("bootstrap", Method, fixed = TRUE)) {

      f0 <- sum(exp(-Lambda) / (1 - exp(-Lambda)))
      Prob <- c(f0, as.numeric(table(y))) / N
      Prob <- cumsum(Prob)
      names(Prob) <- c(0,as.numeric(names(table(y))))
      StrappedStatistic <- NULL

      if (length(y) < 50) {
        Bootnumber <- 10000
      } else if (length(y) < 100) {
        Bootnumber <- 7500
      } else if (length(y) < 500) {
        Bootnumber <- 5000
      } else if (length(y) < 1000) {
        Bootnumber <- 2500
      } else if (length(y) < 2000) {
        Bootnumber <- 2000
      } else if (length(y) < 10000) {
        Bootnumber <- 500
      } else {
        Bootnumber <- 100
      }

      for (z in 1:Bootnumber) {
        U <- stats::runif(N)
        Strap <- NULL

        for (k in names(Prob)) {
          Strap <- c(Strap, length(which(U <= as.numeric(Prob[k]))))
        }
        names(Strap) <- names(Prob)

        m <- length(Strap)

        while (m != 1) {
          Strap[m] <- Strap[m]-Strap[m - 1]
          m <- m - 1
        }

        full <- NULL

        for (k in names(Strap)) {
          full <- c(full, rep(k, Strap[k]))
        }

        full <- full[full != 0]
        Theta <- IRLS(Dependent = as.numeric(full),
                      Covariates = matrix(rep(1, length(full), ncol = 1)),
                      family = Zero_Truncated_Poisson(),
                      start = .1)$Coefficients
        Theta <- exp(matrix(rep(1, length(full), ncol = 1)) %*% Theta)

        StrappedStatistic <- c(StrappedStatistic, sum(1 / (1 - exp(-Theta))))

      }
      N <- mean(StrappedStatistic)
      if (Method == "bootstrapSD") {
        Variation <- stats::var(StrappedStatistic)
        Confidence_Interval <- c(Lower_Bound = max(N - 1.96 * sqrt(Variation), length(y)),
                                 Upper_Bound = N + 1.96 * sqrt(Variation))
      } else if (Method == "bootstrapPerc") {
        Variation <- stats::var(StrappedStatistic)
        Confidence_Interval <- stats::quantile(StrappedStatistic, c(0.025, 0.975))
        names(Confidence_Interval) <- c("Lower_Bound", "Upper_Bound")
      }
      rm(StrappedStatistic);
      rm(Theta);
      rm(full);
      rm(Strap);
      rm(U);
      rm(m);
    }
  } else if (family$family == "ZTNB") {
    alpha <- Overdispersion
    z <- 1 / alpha
    Lambda <- Parameter
    Pr <- 1 - (1 + z * Lambda) ** (- 1 / z)
    N <- sum(1 / Pr)

    if (Method == "analytic") {
    S <- 1 / (1 + z * Lambda)

    Inform <- as.matrix(-Hess(beta))

    BigTheta1 <- sum(((1 + z * Lambda) ** (z - 1) *
    ((1 + z * Lambda) * log(1 + z * Lambda) - z * Lambda)) / (z ** 2))
    BigTheta2 <- -(as.numeric(Lambda * (S ** (1 - 1 / z)) /
    ((1 - (1 / S) ** z) ** 2))) %*% as.matrix(X)
    BigTheta <- matrix(c(BigTheta1, BigTheta2), ncol = 1)

    f1 <-  t(BigTheta) %*% solve(Inform) %*% BigTheta
    f2 <- sum((1 - Pr) / (Pr ** 2))
    Variation <- f1 + f2

    Confidence_Interval <- c(Lower_Bound = max(N - 1.96 * sqrt(Variation), length(y)),
                             Upper_Bound = N + 1.96 * sqrt(Variation))
    names(Confidence_Interval) <- c("Lower_Bound", "Upper_Bound")
    }
  }

  list(Point_estimate = N,
       Variance = Variation,
       Confidence_Interval = Confidence_Interval)
}
