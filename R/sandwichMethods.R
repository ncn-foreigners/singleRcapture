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
#' @return a 
#' @seealso [sandwich::bread()]
#' @importFrom sandwich bread
#' @method bread singleR
#' @exportS3Method
bread.singleR <- function(object,...) {
  return(stats::vcov(object, ...) * as.vector(object$df.residual + length(object$coefficients)))
}