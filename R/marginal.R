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
  y <- table(object$y)[range]
  y <- y[!is.na(y)]
  trcount <- object$trcount
  if(!is.null(onecount)) {
    trcount <- onecount
  }
  
  probFun <- ifelse(object$model$family %in% c("ztpoisson",
                                               "zotpoisson",
                                               "chao",
                                               "zelterman"),
                    function(x, lambda, disp) {stats::dpois(x = x, 
                                                            lambda = lambda)},
                    ifelse(object$model$family %in% c("ztnegbin",
                                                      "zotnegbin"),
                           function(x, lambda, disp) {stats::dnbinom(x = x,
                                                                    mu = lambda,
                                                                    size = exp(-disp))},
                           ifelse(object$model$family %in% c("ztgeom",
                                                             "zotgeom"),
                                  function(x, lambda, disp) {stats::dgeom(x = x,
                                                                          prob = (1 / (1 + lambda)))},
                                  "")))
  ifelse(object$model$family %in% c("zotpoisson",
                                    "zotnegbin",
                                    "zotgeom"),
         probFun2 <- function(x, lambda, disp) {
           (probFun(x = x, lambda = lambda, disp = disp) / 
              (1 - probFun(x = 0, lambda = lambda, disp = disp) - 
               probFun(x = 1, lambda = lambda, disp = disp)))
         },
         probFun2 <- function(x, lambda, disp) {
           (probFun(x = x, lambda = lambda, disp = disp) / 
              (1 - probFun(x = 0, lambda = lambda, disp = disp)))
         }
  )
  
  res <- colSums(t(sapply(object$fitt.values,
          FUN = function(y) {probFun2(x = range,
                                      lambda = y,
                                      disp = object$dispersion)})))
  names(res) <- as.character(range)
  
  if(includeones & (object$model$family %in% c("zotpoisson",
                                               "zotnegbin",
                                               "zotgeom"))) {
    res <- c(trcount, res)
    names(res)[1] <- "1"
    if (!is.null(onecount)) {res[1] <- onecount}
  }
  
  if(includezeros) {
    res <- c(object$populationSize$pointEstimate - length(object$y), res)
    names(res)[1] <- "0"
  }
  res <- structure(list(table = res, y = y, 
                        df = length(y) - 1 - length(object$coefficients)), 
                   class = c("singleRmargin"))
  res
}