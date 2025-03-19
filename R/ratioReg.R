#' @title Estimate population size using ratio regression
#' @author Cyprian Jurkowski, Piotr Chlebicki
#'
#' @description
#' Estimates population size using ratio regression, an alternative method 
#' based on observed capture ratios. This method assumes an exponential 
#' decay pattern in capture probabilities to infer missing individuals.
#'
#' @details
#' **How the function works:**
#' 
#' 1. The function takes a dataset where the `observed` column represents 
#'    how many times an individual was captured.
#' 2. It creates a frequency table (`counts`) to determine how often 
#'    each count appears in the dataset.
#' 3. Using these frequencies, it calculates **capture ratios**:
#'    \deqn{r_k = \frac{f_{k+1}}{f_k}}
#'    where \( f_k \) is the frequency of individuals captured exactly \( k \) times.
#' 4. The function performs a **linear regression** on the logarithm of these ratios:
#'    \deqn{\log(r_k) = \beta_0 + \beta_1 k}
#' 5. The estimated **capture probability** from this regression is then used 
#'    to estimate the **total population size** using:
#'    \deqn{N = \frac{N_{obs}}{1 - P(Y = 0)}}
#' 6. Confidence intervals are computed based on the **variance** of the estimated missing cases.
#'
#' **Key assumptions:**
#' - Capture probabilities follow an exponential decay pattern.
#' - The number of missing individuals can be inferred from observed capture ratios.
#'
#' @param data A data frame containing capture counts (e.g., number of times an individual was observed).
#' @param formula A formula specifying the model structure (e.g., `observed ~ 1`).
#' @param weights Optional prior weights.
#' @param useWeights Logical; whether to apply optional weighting (TODO: Future feature).
#' @param ... Additional arguments.
#' 
#' TODO: Handle optional weighting properly (not yet implemented).
#' TODO: Investigate warnings regarding 'newdata' in predict().
#' TODO: Fix multiple confidence intervals.
#' 
#' @return A list containing:
#' \itemize{
#'   \item{\code{estimates} -- Estimated total population size.}
#'   \item{\code{confidenceInterval} -- Confidence interval for the estimates.}
#'   \item{\code{ratioModel} -- Linear model used for ratio estimation.}
#'   \item{\code{formula} -- The formula used in the regression.}
#'   \item{\code{call} -- The function call.}
#' }
#'
#' @examples
#' # Example 1: Small dataset
#' test_data <- data.frame(observed = c(1, 1, 2, 2, 2, 3, 3, 3, 4, 5))
#' result <- ratioReg(data = test_data, formula = observed ~ 1)
#' print(result$estimates)
#' print(result$confidenceInterval)
#' summary(result$ratioModel)
#'
#' # Example 2: Larger dataset
#' set.seed(123)
#' large_data <- data.frame(observed = sample(1:10, 100, replace = TRUE))
#' result_large <- ratioReg(data = large_data, formula = observed ~ 1)
#' print(result_large$estimates)
#' 
#' @export

ratioReg <- function(data, formula, weights = NULL, useWeights = TRUE, ...) {
  
  modelFrame <- stats::model.frame(formula, data, ...)
  observed <- model.response(modelFrame) |> as.vector()
  
  # Check for negative values
  if (any(observed < 0)) {
    stop("Negative values in 'observed' data are not allowed.")
  }
  
  counts <- table(observed)
  
  # Check for zero counts
  if (any(counts == 0)) {
    stop("Zero counts in the data are not allowed.")
  }
  
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
  
  ff <- log(r) ~ I(seq_len(length(r)))
  model_data <- data.frame(r = r, x = seq_len(length(r))) 
  
  linearModel <- lm(ff, data = model_data, ...)
  
  # TODO: fix problems with structure
  
  new_data <- data.frame(x = seq_len(length(r)))  
  fitRatio <- exp(predict(linearModel, newdata = new_data))  
  
  # TODO: implement weighting option
  
  N <- sum(counts) / (1 - 1 / sum(c(1, cumprod(fitRatio))))
  N <- c(N, sum(counts) + counts["1"] / fitRatio[1])
  names(N) <- c("ht", "reg")
  
  f0 <- N - sum(counts)
  
  variation <- model.matrix(ff, model.frame(ff, new_data))
  variation <- f0^2 * as.vector(variation %*% vcov(linearModel) %*% t(variation))
  
  variation <- variation + counts["1"] * exp(-predict(linearModel, newdata = new_data)) ^ 2
  
  sd <- sqrt(variation)
  sc <- qnorm(p = 1 - .05 / 2)
  
  # TODO: fix computing multiple interval
  
  confidenceInterval <- data.frame(
    lowerBound = pmax(N - sc * sd, sum(counts)),
    upperBound = N + sc * sd
  )
  
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

