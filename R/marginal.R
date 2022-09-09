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
  probFun <- object$model$densityFunction
  
  res <- sapply(1:nrow(object$linear.predictors), 
                FUN = function(x) {object$model$densityFunction(x = range, eta = matrix(object$linear.predictors[x, ], ncol = object$model$parNum))})
  res <- rowSums(res)
  names(res) <- as.character(range)
  
  if(isTRUE(includeones) & (object$model$family %in% c("zotpoisson",
                                                       "zotnegbin",
                                                       "zotgeom"))) {
    res <- c(trcount, res)
    names(res)[1] <- "1"
    y <- c(trcount, y)
    names(y)[1] <- "1"
    if (!is.null(onecount)) {res[1] <- onecount; y[1] <- onecount}
  }
  
  if(isTRUE(includezeros)) {
    res <- c(object$populationSize$pointEstimate - length(object$y), res)
    names(res)[1] <- "0"
  }
  res <- structure(list(table = res, y = y, 
                        df = length(y) - length(object$coefficients) - 1,
                        name = object$model$family), 
                   class = c("singleRmargin"))
  res
}
