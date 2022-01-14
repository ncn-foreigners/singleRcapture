#' Weighted_sd
#'
#' @param x a vector of values from which standard deviation is to be calculated
#' @param w a vector of weights
#' @param ... arguments to be passed from/to other methods
#'
#' @return numeric vector of weighted standard deviation
#' @export
Weighted_Sd <- function(x, w,...){
 Z <- stats::weighted.mean(x, w,...)
 M <- length(w[w != 0])
 sqrt(sum(w * ((x - Z) ** 2)) / ((1-1 / M) * sum(w)))
}
