#' Marginal Frequencies
#'
#' @param object object of class singleR
#' @param includeones Boolean value indicating whether to include estimated number of zero counts
#' @param includezeros Boolean value indicating whether to include one counts in zero one truncated models
#' @param onecount a numeric value indicating number of one counts if null trcount from object will be assumed to be a number one counts
#' @param range optional argument specifying range of selected Y values
#'
#' @return A table with marginal frequencies predicted by the model
#' @export
marginalFreq <- function(object,
                         includeones = TRUE, # matters only for one truncated models
                         includezeros = TRUE,
                         onecount = NULL,
                         range) {
  if (missing(range)) {range <- (min(object$y):max(object$y))}
  y <- table(object$y)[names(table(object$y)) %in% as.character(range)]
  y <- y[!is.na(y)]
  trcount <- object$trcount
  if(!is.null(onecount)) {
    trcount <- onecount
  }
  
  # PMF for truncated distributions:
  probFun <- switch(object$model$family,
  "ztpoisson"  = function(x, lambda, disp) {stats::dpois(x = x, lambda = lambda) / (1 - stats::dpois(x = 0, lambda = lambda))},
  "chao"       = function(x, lambda, disp) {stats::dpois(x = x, lambda = lambda) / (1 - stats::dpois(x = 0, lambda = lambda))},
  "zelterman"  = function(x, lambda, disp) {stats::dpois(x = x, lambda = lambda) / (1 - stats::dpois(x = 0, lambda = lambda))},
  "zotpoisson" = function(x, lambda, disp) {stats::dpois(x = x, lambda = lambda) / (1 - stats::dpois(x = 0, lambda = lambda) - stats::dpois(x = 1, lambda = lambda))},
  "ztnegbin"   = function(x, lambda, disp) {stats::dnbinom(x = x, mu = lambda, size = exp(-disp)) / (1 - stats::dnbinom(x = 0, mu = lambda, size = exp(-disp)))},
  "zotnegbin"  = function(x, lambda, disp) {stats::dnbinom(x = x, mu = lambda, size = exp(-disp)) / (1 - stats::dnbinom(x = 0, mu = lambda, size = exp(-disp)) - stats::dnbinom(x = 1, mu = lambda, size = exp(-disp)))},
  "ztgeom"     = function(x, lambda, disp) {stats::dgeom(x = x, prob = (1 / (1 + lambda))) / (1 - stats::dgeom(x = 0, prob = (1 / (1 + lambda))))},
  "zotgeom"    = function(x, lambda, disp) {stats::dgeom(x = x, prob = (1 / (1 + lambda))) / (1 - stats::dgeom(x = 0, prob = (1 / (1 + lambda))) - stats::dgeom(x = 1, prob = (1 / (1 + lambda))))})

  
  res <- colSums(t(sapply(object$fitt.values$link,
          FUN = function(y) {probFun(x = range,
                                     lambda = y,
                                     disp = object$dispersion)})))
  names(res) <- as.character(range)
  
  if(isTRUE(includeones) & (object$model$family %in% c("zotpoisson",
                                                       "zotnegbin",
                                                       "zotgeom"))) {
    res <- c(trcount, res)
    names(res)[1] <- "1"
    if (!is.null(onecount)) {res[1] <- onecount}
  }
  
  if(isTRUE(includezeros)) {
    res <- c(object$populationSize$pointEstimate - length(object$y), res)
    names(res)[1] <- "0"
  }
  res <- structure(list(table = res, y = y, 
                        df = length(y) - length(object$coefficients),
                        name = object$model$family), 
                   class = c("singleRmargin"))
  res
}
