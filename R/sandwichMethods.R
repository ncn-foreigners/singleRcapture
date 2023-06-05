# ALL CODE IN THIS FILE WAS DEVELOPED BASED ON CODE IN sandwich PACKAGE

# Same story as with bread
#' @importFrom sandwich estfun
#' @method estfun singleR
#' @rdname vcovHC.singleR
#' @exportS3Method
estfun.singleR <- function(x,...) {
  Y <- if (is.null(x$y)) stats::model.response(model.frame(x)) else x$y
  Y <- Y[x$which$reg]
  X <- stats::model.matrix(x, type = "vlm")
  beta <- stats::coef(x)
  wts <- stats::weights(x)
  if (is.null(wts)) wts <- rep(1, length(Y))
  res <- x$model$makeMinusLogLike(y = Y, X = X, NbyK = TRUE, deriv = 1)(beta)
  colnames(res) <- names(beta)
  rownames(res) <- rownames(X)
  res
}

# this literally does the same as every other bread method, there's no need for
# separate documentation just link it to vcovHC
#' @importFrom sandwich bread
#' @method bread singleR
#' @rdname vcovHC.singleR
#' @exportS3Method
bread.singleR <- function(x,...) {
  stats::vcov(x, ...) * as.vector(x$dfResidual + length(coef(x)))
}

#' @title Heteroscedasticity-Consistent Covariance Matrix Estimation for singleR class
#' @author Piotr Chlebicki, Maciej Beręsewicz
#' 
#' @description S3 method for \code{vcovHC} to handle \code{singleR} class objects. 
#' Works exactly like \code{vcov.default} the only difference being that this method handles vector generalised linear models.
#' Updating the covariance matrix in variance/standard error estimation for population size estimator can be done via [singleRcapture::redoPopEstimation()]
#'
#' @param x a fitted \code{singleR} class object.
#' @param type a character string specifying the estimation type, same as in \code{sandwich::vcovHC.default}. HC3 is the default value.
#' @param omega a vector or a function depending on the arguments residuals (i.e. the derivative of log-likelihood with respect to each linear predictor), diaghat (the diagonal of the corresponding hat matrix) and df (the residual degrees of freedom), same as in \code{sandwich::vcovHC.default}.
#' @param sandwich logical. Should the sandwich estimator be computed? If set to FALSE only the meat matrix is returned. Same as in [sandwich::vcovHC()]
#' @param ... for \code{vcovHC} additional optional arguments passed to the following functions:
#' \itemize{
#'   \item \code{estfun} -- for empirical estimating functions.
#'   \item \code{hatvalues} -- for diagonal elements of projection matrix.
#'   \item \code{sandwich} -- only if \code{sandwich} argument in function call was set to \code{TRUE}.
#'   \item \code{vcov} -- when calling \code{bread} internally.
#' }
#'
#' @return Variance-covariance matrix estimation corrected for heteroscedasticity of regression errors.
#' @seealso [sandwich::vcovHC()] [singleRcapture::redoPopEstimation()]
#' @examples 
#' set.seed(1)
#' N <- 10000
#' gender <- rbinom(N, 1, 0.2)
#' eta <- -1 + 0.5*gender
#' counts <- rpois(N, lambda = exp(eta))
#' df <- data.frame(gender, eta, counts)
#' df2 <- subset(df, counts > 0)
#' mod1 <-  estimatePopsize(
#'   formula = counts ~ 1 + gender, 
#'   data = df2, 
#'   model = "ztpoisson", 
#'   method = "optim", 
#'   popVar = "analytic"
#' )
#' require(sandwich)
#' HC <- sandwich::vcovHC(mod1, type = "HC4")
#' Fisher <- vcov(mod1, "Fisher") # variance covariance matrix obtained from 
#' #Fisher (expected) information matrix
#' HC
#' Fisher
#' # usual results
#' summary(mod1)
#' # updated results
#' summary(mod1, cov = HC,
#' popSizeEst = redoPopEstimation(mod1, cov = HC))
#' # estimating equations
#' mod1_sims <- sandwich::estfun(mod1)
#' head(mod1_sims)
#' # bread method
#' all(vcov(mod1, "Fisher") * nrow(df2) == sandwich::bread(mod1, type = "Fisher"))
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
  beta <- x$coefficients
  X <- model.matrix(x, "vlm")
  n <- nrow(X)
  k <- ncol(X)
  df <- n - k
  hat <- as.vector(hatvalues(x, ...))
  Y <- if (is.null(x$y)) stats::model.response(model.frame(x)) else x$y
  Y <- Y[x$which$reg] # only choose units which appear in regression
  res <- as.vector(x$model$makeMinusLogLike(y = Y, X = X, vectorDer = TRUE, der = 1)(beta))
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
        delta <- pmin(n * diaghat/p, pmax(4, n * k * max(diaghat)/p))
        residuals^2 / sqrt((1 - diaghat)^delta)
      }
    })
  }
  if (is.function(omega)) 
    omega <- omega(res, hat, df)
  rval <- sqrt(omega) * X
  rval <- crossprod(rval)/n
  if (sandwich) 
    rval <- sandwich(x, meat. = rval, ...)
  rval
}