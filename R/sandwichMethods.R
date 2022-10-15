# TODO::
## - boot (bs)
## - HC in CL
#' Empirical Estimating Functions
#' 
#' An S3method for \code{sandwich::estfun} to handle \code{singleR} objects.
#' 
#' @param object an object representing a fitted model.
#' @param ... additional optional arguments.
#' @return A \code{matrix} with \code{n} rows and \code{k} columns where \code{k} denotes number of variables.
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

#' Extracting bread matrix for singleR class
#' 
#' An S3class for \code{sandwich::bread} to handle \code{singleR} objects. This function was developed based on \code{sandwich:::bread.glm}
#' 
#' @param object an object representing a fitted model.
#' @param ... additional optional arguments passed to the following functions:
#' \itemize{
#'   \item \code{stats::vcov} -- for extracting "the usual" variance-covariance matrix, \code{vcov.singleR} has one additional argument \code{type} with values \code{"Fisher"} \code{"observedInform"}, defaults to the one specified in \code{control.pop.var} specified in call for object.
#' }
#' @return A bread matrix, i.e. a hessian based estimation of variance-covariance matrix scaled by degrees of freedom.
#' @seealso [sandwich::bread()] [singleRcapture::control.pop.var()]
#' @examples 
#' Model <- estimate_popsize(
#' formula = capture ~ ., 
#' data = netherlandsimmigrant, 
#' model = ztpoisson, 
#' method = "robust")
#' sandwich::bread(Model)
#' vcov(Model)
#' # This function just scales.
#' all(vcov(Model) * nrow(netherlandsimmigrant) == sandwich::bread(Model))
#' # We can choose Fisher information matrix instead of default observed information matrix.
#' vcov(Model, "Fisher")
#' sandwich::bread(Model, type = "Fisher")
#' all(vcov(Model, "Fisher") * nrow(netherlandsimmigrant) == sandwich::bread(Model, type = "Fisher"))
#' @importFrom sandwich bread
#' @method bread singleR
#' @exportS3Method
bread.singleR <- function(object,...) {
  return(stats::vcov(object, ...) * as.vector(object$df.residual + length(object$coefficients)))
}

#' Heteroscedasticity-Consistent Covariance Matrix Estimation for singleR class
#' 
#' @description S3 method for \code{vcovHC} to handle \code{singleR} class objects. 
#' Works exactly like \code{vcov.default} the only difference being that this method handles vector generalised linear models.
#' Updating the covariance matrix in variance/standard error estimation for population size estimator can be done via [singleRcapture::redoPopEstimation()]
#'
#' @param x a fitted \code{singleR} class object.
#' @param type a character string specifying the estimation type, same as in \code{sandwich::vcovHC.default}. HC3 is the default value.
#' @param omega a vector or a function depending on the arguments residuals (i.e. the derivative of log-likelihood with respect to each linear predictor), diaghat (the diagonal of the corresponding hat matrix) and df (the residual degrees of freedom), same as in \code{sandwich::vcovHC.default}.
#' @param sandwich logical. Should the sandwich estimator be computed? If set to FALSE only the meat matrix is returned. Same as in [sandwich::vcovHC()]
#' @param ... additional optional arguments passed to the following functions:
#' \itemize{
#'   \item estfun -- for empirical estimating functions.
#'   \item hatvalues -- for diagonal elements of projection matrix.
#'   \item sandwich -- only if \code{sandwich} argument in function call was set to \code{TRUE}.
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
#' mod1 <-  estimate_popsize(formula = counts ~ 1 + gender, data = df2, 
#' model = "ztpoisson", method = "mle", pop.var = "analytic")
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
  hat <- as.vector(hatvalues(x, ...))
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
  return(rval)
}
