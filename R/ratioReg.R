#' @title Estimate population size using ratio regression
#' @author Cyprian Jurkowski, Piotr Chlebicki
#' @description Estimates population size using ratio regression, an alternative
#' method based on ratios of observed capture counts.
#' 
#' ## TODO: Handle optional weighting properly (not yet implemented).
#' ## TODO: Investigate warnings regarding 'newdata' in predict().
#' ## TODO: Fix multiple confidence intervals.
#' 
#' @param data A data frame containing capture counts.
#' @param formula A formula specifying the model structure.
#' @param weights Optional prior weights. # TODO: Implement this in calculations.
#' @param useWeights Logical; whether to apply optional weighting.
#' @param ... Additional arguments.
#' @return A list containing population size estimates and a confidence interval.
#' @export

ratioReg <- function(data, formula, weights = NULL, useWeights = TRUE, ...) {
  
  # TODO: Debug and investigate warnings regarding 'newdata' in predict()
  
  modelFrame <- stats::model.frame(formula, data, ...)
  observed <- model.response(modelFrame) |> as.vector()
  
  counts <- table(observed)
  
  if (length(counts) < 2) {
    stop("Ratio regression requires at least two observed capture counts.")
  }
  
  r <- sapply(1:(max(observed)-1), function(x) {
    if (as.character(x) %in% names(counts) && as.character(x+1) %in% names(counts)) {
      return(counts[as.character(x+1)] / counts[as.character(x)])
    } else {
      return(NA)
    }
  }) |> as.vector()
  
  r <- na.omit(r)
  
  # TODO: Investigate cases where r becomes empty or too short for fitting a model
  
  ff <- log(r) ~ I(seq_len(length(r)))
  model_data <- data.frame(r = r, x = seq_len(length(r))) 
  linearModel <- lm(ff, data = model_data, ...)
  
  new_data <- data.frame(x = seq_len(length(r)))
  fitRatio <- exp(predict(linearModel, newdata = new_data))
  
  N <- sum(counts) / (1 - 1 / sum(c(1, cumprod(fitRatio))))
  N <- c(N, sum(counts) + counts["1"] / fitRatio[1])
  names(N) <- c("ht", "reg")
  
  f0 <- N - sum(counts)
  
  variation <- model.matrix(ff, model.frame(ff, new_data))
  variation <- f0^2 * as.vector(variation %*% vcov(linearModel) %*% t(variation))
  
  # TODO: Investigate variance calculation for small sample sizes
  
  variation <- variation + counts["1"] * exp(-predict(linearModel, newdata = new_data)) ^ 2
  
  sd <- sqrt(variation)
  sc <- qnorm(p = 1 - .05 / 2)
  
  # TODO: Fix multiple confidence intervals.
  
  confidenceInterval <- data.frame(
    lowerBound = pmax(N - sc * sd, sum(counts)),
    upperBound = N + sc * sd
  )
  
  # TODO: Handle optional weighting in calculations
  
  structure(
    list(
      estimates = N,
      confidenceInterval = confidenceInterval,
      ratioModel = linearModel,
      formula = formula,
      call = match.call()
    ),
    class = c("singleRRatioReg", "singleR")
  )
}
