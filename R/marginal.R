#' @title Observed and fitted marginal Frequencies
#' @author Piotr Chlebicki
#' \loadmathjax
#' @description A function that given a fitted \code{singleR} class object 
#' computed marginal frequencies by as sum of probability density functions
#' for each unit in data at each point i.e. kth element of marginal frequency
#' table is given by \mjeqn{\sum_{j=1}^{N_{obs}}\mathbb{P}(Y_{j}=k|
#' \eta_{j})}{sum_j=1^N_obs P(Y_j=k|eta_{j})}. For k=0 only 
#' (if specified at call) they are computed as \mjeqn{\hat{N}-N_{obs}}{N-N_obs} 
#' because \mjeqn{\boldsymbol{f}_{0}}{f_0} is assumed to the unobserved part of the 
#' studied population.
#' 
#' These frequencies are useful in diagnostics for count data regression, such
#' as assessment of fit.
#'
#' @param object Object of \code{singleR} class.
#' @param includeones Logical value indicating whether to include the estimated number of zero counts.
#' @param includezeros Logical value indicating whether to include one counts in the zero-one truncated models.
#' @param onecount A numeric value indicating number of one counts if null \code{trcount} from object will be assumed to be a number one counts.
#' @param range Optional argument specifying range of selected Y values.
#'
#' @return A list with observed name of the fitted model family degrees of freedom and observed and fitted marginal frequencies.
#' @export
marginalFreq <- function(object,
                         includeones = TRUE, # matters only for zero one truncated models
                         includezeros = TRUE,
                         onecount = NULL,
                         range) {
  y <- if (is.null(object$y)) stats::model.response(model.frame(object)) else object$y
  if (missing(range)) {range <- (min(y):max(y))}
  y <- table(y)[names(table(y)) %in% as.character(range)]
  y <- y[!is.na(y)]
  trcount <- object$trcount
  if(!is.null(onecount)) {
    trcount <- onecount
  }
  
  # PMF for truncated distributions:
  probFun <- object$model$densityFunction
  
  res <- sapply(
    1:nrow(object$linear.predictors), 
    FUN = function(x) {object$model$densityFunction(
      x = range, 
      eta = matrix(object$linear.predictors[x, ], 
                   ncol = object$model$parNum)
    )}
  )
  res <- rowSums(res)
  names(res) <- as.character(range)
  
  if(isTRUE(includeones) & grepl(object$model$family, pattern = "^zot")) {
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
  res <- structure(list(
    table = res, y = y, 
    df = length(y) - length(object$coefficients) - 1,
    name = object$model$family
  ),
  class = c("singleRmargin"))
  res
}
