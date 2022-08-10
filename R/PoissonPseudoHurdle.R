#' Pseudohurdle poisson
#' 
#' @return A object of class "family" containing objects \cr
#' makeMinusLogLike(y,X) - for creating negative likelihood function \cr
#' makeGradient(y,X) - for creating gradient function \cr
#' makeHessian(X) - for creating hessian \cr
#' linkfun - a link function to connect between linear predictor and model parameter in regression and a name of link function\cr
#' linkinv - an inverse function of link \cr
#' Dlink - a 1st derivative of link function \cr
#' mu.eta,Variance - Expected Value and Variance \cr
#' valedmu, valideta - for checking if regression arguments and valid\cr
#' family - family name\cr
#' Where: \cr
#' y is a vector of observed values \cr
#' X is a matrix / data frame of covariates
#' @importFrom lamW lambertW0
#' @export
# pshurdlepoisson <- function() {
#   link <- log
#   invlink <- exp
#   dlink <- function(lambda) {
#     1 / lambda
#   }
#   
#   mu.eta <- function(eta, type = "trunc", theta, ...) {
#     #TODO
#     lambda <- invlink(eta)
#     omega <- theta / (1 + theta)
#     switch (type,
#             "nontrunc" = (omega + (1 - omega) * lambda / (exp(lambda) - 1) + (1 - omega) * lambda) * (1 - exp(-lambda)),
#             "trunc" = omega + (1 - omega) * lambda / (exp(lambda) - 1) + (1 - omega) * lambda
#     )
#   }
#   
#   variance <- function(mu, type = "nontrunc", theta, ...) {
#     #TODO
#     omega <- theta / (1 + theta)
#     switch (type,
#             "nontrunc" = ((omega + (1 - omega) * mu / (exp(mu) - 1) + (1 - omega) * mu + (mu ** 2) * (1 - omega) / (1 - exp(-mu))) * (1 - exp(-mu))- 
#                             ((omega + (1 - omega) * mu / (exp(mu) - 1) + (1 - omega) * mu) * (1 - exp(-mu))) ** 2),
#             "trunc" = ((omega + (1 - omega) * mu / (exp(mu) - 1) + (1 - omega) * mu + (mu ** 2) * (1 - omega) / (1 - exp(-mu))) - 
#                          ((omega + (1 - omega) * mu / (exp(mu) - 1) + (1 - omega) * mu) ** 2))
#     )
#   }
#   
#   Wfun <- function(prior, eta, theta, ...) {
#     #TODO
#     lambda <- exp(eta)
#     z <- 1 - exp(-lambda) # expected for I's
#     (z * (-1 / ((theta + lambda / (exp(lambda) - 1)) ** 2) *
#             (((exp(lambda) - 1 - lambda * exp(lambda)) / ((exp(lambda) - 1) ** 2)) ** 2)
#           + (lambda / (theta + lambda / (exp(lambda) - 1)))
#           * (-lambda * exp(lambda) * ((exp(lambda) - 1) ** 2) 
#              -2 * (exp(lambda) - 1) * exp(lambda) * (exp(lambda) - 1 - lambda * exp(lambda)))
#           / ((exp(lambda) - 1) ** 4)) - (1 - z) * lambda * ((exp(-lambda) + lambda * exp(- lambda) - 1) / ((1 - exp(-lambda)) ** 2)))
#   }
#   
#   funcZ <- function(eta, weight, y, mu, theta, ...) {
#     #TODO
#     lambda <- exp(eta)
#     z <- ifelse(y == 1, y, 0)
#     eta + (z * (1 / (theta + lambda / (exp(lambda) - 1)) *
#                   (exp(lambda) - 1 - lambda * exp(lambda)) / 
#                   ((exp(lambda) - 1) ** 2)) + 
#              ((1 - z) * (y - mu))) / weight#TODO
#   }
#   
#   minusLogLike <- function(y, X, weight = 1, type = "abc") {
#     #TODO
#     y <- as.numeric(y)
#     if (is.null(weight)) {
#       weight <- 1
#     }
#     z <- ifelse(y == 1, y, 0)
#     function(arg) {
#       beta <- arg
#       eta <- as.matrix(X) %*% beta
#       lambda <- exp(eta)
#       -sum(weight * ((1 - prob) * (1 - z) * (y * eta - lambda - log(factorial(y)) -
#                                             log(1 - exp(-lambda) - lambda * exp(-lambda)))))
#     }
#   }
#   
#   gradient <- function(y, X, weight = 1) {
#     #TODO
#     y <- as.numeric(y)
#     if (is.null(weight)) {
#       weight <- 1
#     }
#     z <- ifelse(y == 1, y, 0)
#     function(arg) {
#       beta <- arg[-1]
#       theta <- arg[1]
#       eta <- as.matrix(X) %*% beta
#       lambda <- exp(eta)
#       mu <- lambda / (1 - exp(-lambda))
#       c(sum(weight * (- 1 / (1 + theta) + z / (theta + lambda / (exp(lambda) - 1)))), # theta derivative
#         t(as.matrix(X)) %*% ((z * (1 / (theta + lambda / (exp(lambda) - 1)) *  # beta derivative
#                                      (exp(lambda) - 1 - lambda * exp(lambda)) / 
#                                      ((exp(lambda) - 1) ** 2)) + 
#                                 ((1 - z) * weight * (y - mu))) * weight))
#     }
#   }
#   
#   hessian <- function(y, X, weight = 1) {
#     #TODO
#     if (is.null(weight)) {
#       weight <- 1
#     }
#     z <- ifelse(y == 1, y, 0)
#     function(arg) {
#       beta <- arg[-1]
#       theta <- arg[1]
#       res <- matrix(1, nrow = length(arg), ncol = length(arg), dimnames = list(names(arg), names(arg)))
#       lambda <- exp(as.matrix(X) %*% beta)
#       eml <- exp(-lambda)
#       coefficient <- (1 / (1 - eml) - lambda * eml / ((1 - eml) ** 2))
#       
#       dmu <- weight * as.numeric(coefficient)
#       dlam <- as.matrix(X * as.numeric(lambda))
#       
#       # Theta^2 derivative
#       
#       G00 <- sum((1 / ((1 + theta) ** 2) - z / ((theta + lambda / (exp(lambda) - 1)) ** 2)) * weight)
#       
#       # Theta beta derivative
#       
#       G01 <- (t(as.matrix(X)) %*% ((-z * ((theta + lambda / (exp(lambda) - 1)) ** -2)
#                                     * (exp(lambda) - 1 - lambda * exp(lambda))
#                                     / ((exp(lambda) - 1) ** 2)) * weight))
#       
#       # Beta^2 derivative
#       G11 <- ((t(as.matrix(X) *
#                    as.numeric(weight * z * (-1 / ((theta + lambda / (exp(lambda) - 1)) ** 2) *
#                                               (((exp(lambda) - 1 - lambda * exp(lambda)) / ((exp(lambda) - 1) ** 2)) ** 2)
#                                             + (lambda / (theta + lambda / (exp(lambda) - 1)))
#                                             * (-lambda * exp(lambda) * ((exp(lambda) - 1) ** 2) 
#                                                -2 * (exp(lambda) - 1) * exp(lambda) * (exp(lambda) - 1 - lambda * exp(lambda)))
#                                             / ((exp(lambda) - 1) ** 4)))) %*% as.matrix(X))
#               -((t(as.matrix(X) * ((1 - z) * dmu))) %*% dlam))
#       
#       
#       res[1, 1] <- G00
#       res[-1, -1] <- G11
#       res[1, -1] <- G01
#       res[-1, 1] <- G01
#       
#       res
#     }
#   }
#   
#   validmu <- function(mu) {
#     #TODO
#     (sum(!is.finite(mu)) == 0) && all(0 < mu)
#   }
#   
#   dev.resids <- function(y, mu, wt, theta, ...) {
#     #TODO
#     NULL
#   }
#   
#   pointEst <- function (pw, lambda, contr = FALSE, theta, ...) {
#     #TODO
#     N <- pw / (1 - exp(-lambda))
#     if(!contr) {
#       N <- sum(N)
#     }
#     N
#   }
#   
#   popVar <- function (pw, lambda, cov, X, theta, ...) {
#     #TODO
#     X <- as.data.frame(X)
#     ml <- (1 - exp(-lambda)) ** 2
#     
#     f1 <- colSums(-X * pw * (exp(log(lambda) - lambda) / ml))
#     f1 <- c(0, t(f1)) %*% as.matrix(cov) %*% c(0, f1)
#     
#     f2 <- sum(pw * exp(-lambda) / ml)
#     
#     f1 + f2
#   }
#   
#   structure(
#     list(
#       makeMinusLogLike = minusLogLike,
#       makeGradient = gradient,
#       makeHessian = hessian,
#       linkfun = link,
#       linkinv = invlink,
#       dlink = dlink,
#       mu.eta = mu.eta,
#       link = "log",
#       valideta = function (eta) {TRUE},
#       variance = variance,
#       Wfun = Wfun,
#       funcZ = funcZ,
#       dev.resids = dev.resids,
#       validmu = validmu,
#       pointEst = pointEst,
#       popVar= popVar,
#       family = "hurdlepoisson",
#       parNum = 2
#     ),
#     class = "family"
#   )
# }
