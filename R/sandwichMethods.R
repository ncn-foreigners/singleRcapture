# TODO::
## - HC
## - boot (bs)
## - HC in CL

#' estfun
#' 
#' An S3class for \code{sandwich::estfun} to handle \code{singleR} objects. This function was developed based on \code{countreg::estfun.zerotrunc}
#' 
#' @param object an object representing a fitted model.
#' @param ... additional optional arguments.
#' @return a \code{matrix} with \code{n} rows and \code{k} columns where \code{k} denotes number of variables.
#' @seealso [sandwich::estfun()]
#' @examples 
#' set.seed(1)
#' N <- 10000
#' gender <- rbinom(N, 1, 0.2)
#' eta <- -1 + 0.5*gender
#' counts <- rpois(N, lambda = exp(eta))
#' df <- data.frame(gender, eta, counts)
#' df2 <- subset(df, counts > 0)
#' mod1 <-  estimate_popsize(formula = counts ~ 1 + gender, data = df2, 
#' model = "ztpoisson", method = "mle", pop.var = "analytic")
#' mod1_sims <- sandwich::estfun(mod1)
#' head(mod1_sims) 
#' @importFrom sandwich estfun
#' @method estfun singleR
#' @exportS3Method
estfun.singleR <- function(object,...) {
  #if (object$model$parNum != 1) stop("For now only models with single linear predictor have estimating functions and by extension sandwitch functionalities implemented.")
  Y <- if (is.null(object$y)) stats::model.response(model.frame(object)) else object$y
  Y <- Y[object$which$reg]
  X <- stats::model.matrix(object, type = "vlm")[[1]]
  beta <- stats::coef(object)
  wts <- stats::weights(object)
  if (is.null(wts)) wts <- rep(1, length(Y))
  res <- object$model$makeGradient(y = Y, X = X, NbyK = TRUE)(beta)
  colnames(res) <- names(beta)
  rownames(res) <- rownames(X)
  res
}

#' bread
#' 
#' An S3class for \code{sandwich::bread} to handle \code{singleR} objects. This function was developed based on \code{sandwich:::bread.glm}
#' 
#' @param object an object representing a fitted model.
#' @param ... additional optional arguments.
#' @return A bread matrix
#' @seealso [sandwich::bread()]
#' @importFrom sandwich bread
#' @method bread singleR
#' @exportS3Method
bread.singleR <- function(object,...) {
  return(stats::vcov(object, ...) * as.vector(object$df.residual + length(object$coefficients)))
}

#' Title
#'
#' @param x TODO
#' @param type TODO
#' @param omega TODO
#' @param sandwich TODO
#' @param ... TODO
#'
#' @return vcov for hc
#' @importFrom sandwich vcovHC
#' @method vcovHC singleR
#' @exportS3Method
vcovHC.singleR <- function(x, 
                           type = c("HC3", "const", "HC", 
                                    "HC0", "HC1", "HC2", 
                                    "HC4", "HC4m", "HC5"), 
                           omega = NULL, 
                           sandwich = TRUE, 
                           ...) {
  type <- match.arg(type)
  estfun <- estfun(x, ...)
  Y <- x$y
  beta <- x$coefficients
  X <- model.matrix(x, "vlm")
  X <- X[[1]]
  n <- nrow(X)
  k <- ncol(X)
  df <- n - k
  hat <- as.vector(hatvalues(x))
  res <- as.vector(x$model$makeGradient(y = Y, X = X, vectorDer = TRUE)(beta))
  if (is.null(omega)) {
    if (type == "HC") 
      type <- "HC0"
    switch(type, const = {
      omega <- function(residuals, diaghat, df) rep(1, length(residuals)) * sum(residuals^2)/df
    }, HC0 = {
      omega <- function(residuals, diaghat, df) residuals^2
    }, HC1 = {
      omega <- function(residuals, diaghat, df) residuals^2 * 
        length(residuals)/df
    }, HC2 = {
      omega <- function(residuals, diaghat, df) residuals^2/(1 - diaghat)
    }, HC3 = {
      omega <- function(residuals, diaghat, df) residuals^2/(1 - diaghat)^2
    }, HC4 = {
      omega <- function(residuals, diaghat, df) {
        n <- length(residuals)
        p <- as.integer(round(sum(diaghat), digits = 0))
        delta <- pmin(4, n * diaghat/p)
        residuals^2/(1 - diaghat)^delta
      }
    }, HC4m = {
      omega <- function(residuals, diaghat, df) {
        gamma <- c(1, 1.5)
        n <- length(residuals)
        p <- as.integer(round(sum(diaghat), digits = 0))
        delta <- pmin(gamma[1], n * diaghat/p) + pmin(gamma[2], 
                                                      n * diaghat/p)
        residuals^2/(1 - diaghat)^delta
      }
    }, HC5 = {
      omega <- function(residuals, diaghat, df) {
        k <- 0.7
        n <- length(residuals)
        p <- as.integer(round(sum(diaghat), digits = 0))
        delta <- pmin(n * diaghat/p, pmax(4, n * k * 
                                            max(diaghat)/p))
        residuals^2/sqrt((1 - diaghat)^delta)
      }
    })
  }
  if (is.function(omega)) 
    omega <- omega(res, hat, df)
  rval <- sqrt(omega) * X
  rval <- crossprod(rval)/n
  if (sandwich) 
    rval <- sandwich(x, meat. = rval, ...)
  return(rval)
}