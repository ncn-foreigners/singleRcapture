#' @title Estimate Population Size Using Ratio Regression
#' @author Cyprian Jurkowski, Piotr Chlebicki, Maciej Beręsewicz
#'
#' @description
#' Fits the ratio-regression approach for single-source capture-recapture
#' count data as a standalone estimator. The regression is performed on
#' neighbouring marginal-frequency ratios rather than on unit-level counts,
#' which is why this method is separate from the `singleRmodels()` family
#' framework used by [estimatePopsize()].
#'
#' @details
#' Let \eqn{f_k} denote the observed frequency of units captured exactly
#' \eqn{k} times. Ratio regression models the neighbouring ratios
#' \eqn{r_k = f_{k + 1} / f_k}. The default model uses
#' \eqn{\log(r_k) = \beta_0 + \beta_1 \log(k + 1)} and weighted least squares
#' with weights \eqn{1 / (1 / f_k + 1 / f_{k + 1})}.
#'
#' Two published models are supported:
#'
#' - `M0`: ratio regression without one-inflation term
#' - `M1`: ratio regression with an additional one-inflation indicator
#'   `I(k == 1)`
#'
#' When `model = "auto"`, both models are fitted and compared by `criterion`.
#' Population size can then be extracted either through the Horvitz-Thompson
#' (`HT`) or semi-parametric (`SM`) estimator, with `estimator` selecting which
#' one is treated as primary by [popSizeEst()], [summary()], and [plot()].
#'
#' Bootstrap intervals are percentile intervals. For fixed `model = "M0"` or
#' `model = "M1"` the bootstrap uses zero-count imputation based on the
#' primary estimator. For `model = "auto"` it uses the single-bootstrap
#' reselection scheme where model selection is repeated inside each bootstrap
#' sample.
#'
#' @param formula a formula identifying the observed count column. It must use
#'   `~ 1`; for example `observed ~ 1`.
#' @param data a data frame or an object coercible to `data.frame`.
#' @param ratioFormula a one-sided formula used for the ratio regression on the
#'   derived ratio-level data. The default is `~ log(k + 1)`.
#' @param model which ratio-regression model to use: `M0`, `M1`, or `auto` to
#'   select between them by `criterion`.
#' @param criterion information criterion used when `model = "auto"`.
#' @param estimator which population-size estimator should be treated as the
#'   primary one.
#' @param weights optional prior weights interpreted as multiplicities of
#'   observed units.
#' @param subset optional logical expression selecting rows from `data`.
#' @param confint either `"bootstrap"` for percentile bootstrap intervals or
#'   `"none"` to skip interval estimation.
#' @param B number of bootstrap replicates when `confint = "bootstrap"`.
#' @param seed optional random seed used for bootstrap sampling.
#' @param maxCount optional upper support bound. Counts larger than `maxCount`
#'   are discarded before fitting.
#' @param ... additional arguments passed to [stats::model.frame()].
#'
#' @return
#' An object of class `singleRRatioReg` containing the fitted ratio-regression
#' models, model-selection criteria, derived ratio-level data, observed
#' frequencies, and population-size estimates.
#'
#' @examples
#' toy_counts <- data.frame(
#'   observed = c(rep(1, 50), rep(2, 30), rep(3, 15), rep(4, 6), rep(5, 2))
#' )
#'
#' toy_fit <- ratioReg(
#'   observed ~ 1,
#'   data = toy_counts,
#'   model = "auto",
#'   estimator = "SM",
#'   confint = "none"
#' )
#'
#' summary(toy_fit)
#' popSizeEst(toy_fit)
#' plot(toy_fit)
#'
#' toy_tab <- aggregate(list(weight = rep(1, nrow(toy_counts))),
#'                      by = list(observed = toy_counts$observed),
#'                      FUN = sum)
#' ratioReg(
#'   observed ~ 1,
#'   data = toy_tab,
#'   weights = toy_tab$weight,
#'   model = "M0",
#'   confint = "none"
#' )
#'
#' @seealso [estimatePopsize()] [popSizeEst()]
#' @export
ratioReg <- function(formula,
                     data,
                     ratioFormula = ~ log(k + 1),
                     model = c("auto", "M0", "M1"),
                     criterion = c("AIC", "BIC"),
                     estimator = c("SM", "HT"),
                     weights = NULL,
                     subset = NULL,
                     confint = c("bootstrap", "none"),
                     B = 1000,
                     seed = NULL,
                     maxCount = NULL,
                     ...) {
  model <- match.arg(model)
  criterion <- match.arg(criterion)
  estimator <- match.arg(estimator)
  confint <- match.arg(confint)

  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 0) {
    stop("Argument B must be a single non-negative number.")
  }
  B <- as.integer(B)
  if (identical(confint, "bootstrap") && B < 1L) {
    stop("Argument B must be at least 1 when confint = 'bootstrap'.")
  }

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))) {
    stop("Argument seed must be NULL or a single finite number.")
  }

  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }

  singleRratioRegInternalCheckPrimaryFormula(formula)
  ratioFormula <- singleRratioRegInternalCheckRatioFormula(ratioFormula)

  fullDataRows <- rownames(data)
  modelFrameArgs <- list(
    formula = formula,
    data = data,
    ...
  )
  if (!is.null(subset)) {
    modelFrameArgs$subset <- substitute(subset)
  }
  modelFrame <- do.call(stats::model.frame, modelFrameArgs)
  observed <- stats::model.response(modelFrame)

  if (NCOL(observed) > 1L) {
    stop("ratioReg() supports only a single response variable.")
  }

  observed <- as.numeric(observed)
  singleRratioRegInternalCheckObserved(observed)

  priorWeights <- singleRratioRegInternalAlignWeights(
    weights = weights,
    modelFrame = modelFrame,
    data = data,
    fullDataRows = fullDataRows
  )

  counts <- singleRratioRegInternalBuildCounts(
    observed = observed,
    weights = priorWeights,
    maxCount = maxCount
  )

  singleRratioRegInternalFitCounts(
    counts = counts,
    formula = formula,
    ratioFormula = ratioFormula,
    call = match.call(),
    requestedModel = model,
    criterion = criterion,
    primaryEstimator = estimator,
    confintType = confint,
    B = B,
    seed = seed
  )
}

singleRratioRegInternalCheckPrimaryFormula <- function(formula) {
  tt <- stats::terms(formula)
  if (attr(tt, "response") == 0L) {
    stop("Argument formula must include a response variable.")
  }

  if (length(attr(tt, "term.labels")) > 0L || attr(tt, "intercept") == 0L) {
    stop(
      "Argument formula must use '~ 1'. ",
      "Use ratioFormula for the ratio regression terms."
    )
  }

  invisible(formula)
}

singleRratioRegInternalCheckRatioFormula <- function(ratioFormula) {
  tt <- stats::terms(ratioFormula)
  if (attr(tt, "response") != 0L) {
    stop("Argument ratioFormula must be a one-sided formula.")
  }
  ratioFormula
}

singleRratioRegInternalCheckObserved <- function(observed) {
  if (any(!is.finite(observed))) {
    stop("Observed counts must be finite.")
  }
  if (any(observed <= 0)) {
    stop("Observed counts must be strictly positive.")
  }
  if (any(abs(observed - round(observed)) > sqrt(.Machine$double.eps))) {
    stop("Observed counts must be whole numbers.")
  }
  invisible(observed)
}

singleRratioRegInternalAlignWeights <- function(weights,
                                                modelFrame,
                                                data,
                                                fullDataRows) {
  if (is.null(weights)) {
    return(rep(1, nrow(modelFrame)))
  }

  weights <- as.numeric(weights)
  modelFrameRows <- rownames(modelFrame)
  dataRows <- rownames(data)

  if (length(weights) == 1L) {
    priorWeights <- rep(weights, nrow(modelFrame))
  } else if (length(weights) == nrow(modelFrame)) {
    priorWeights <- weights
  } else if (length(weights) == nrow(data)) {
    alignedRows <- match(modelFrameRows, dataRows)
    if (anyNA(alignedRows)) {
      stop("Unable to align argument weights with rows used in ratioReg().")
    }
    priorWeights <- weights[alignedRows]
  } else if (length(weights) == length(fullDataRows)) {
    alignedRows <- match(modelFrameRows, fullDataRows)
    if (anyNA(alignedRows)) {
      stop("Unable to align argument weights with rows used in ratioReg().")
    }
    priorWeights <- weights[alignedRows]
  } else {
    stop(
      "Argument weights must have length 1, the number of rows used in fitting, ",
      "the number of rows in data, or the number of rows in the original data."
    )
  }

  if (any(!is.finite(priorWeights)) || any(priorWeights < 0)) {
    stop("Argument weights must be finite and non-negative.")
  }

  priorWeights
}

singleRratioRegInternalBuildCounts <- function(observed, weights, maxCount = NULL) {
  counts <- stats::aggregate(
    x = list(weight = weights),
    by = list(observed = observed),
    FUN = sum
  )
  counts <- counts[order(counts$observed), , drop = FALSE]

  if (!is.null(maxCount)) {
    if (length(maxCount) != 1L || !is.finite(maxCount)) {
      stop("Argument maxCount must be a single finite number.")
    }
    maxCount <- as.integer(maxCount)
    if (maxCount < 2L) {
      stop("Argument maxCount must be at least 2.")
    }
    counts <- subset(counts, observed <= maxCount)
  }

  if (!nrow(counts)) {
    stop("No positive counts remained after applying maxCount.")
  }

  supportMax <- max(counts$observed)
  fullCounts <- rep(0, supportMax)
  names(fullCounts) <- as.character(seq_len(supportMax))
  fullCounts[as.character(counts$observed)] <- counts$weight
  if (sum(fullCounts) <= 0) {
    stop("Weighted observed frequencies must sum to a positive value.")
  }
  fullCounts
}

singleRratioRegInternalBuildRatioData <- function(counts) {
  if (length(counts) < 2L) {
    stop("ratioReg() requires support up to at least count 2.")
  }

  k <- seq_len(length(counts) - 1L)
  keep <- counts[k] > 0 & counts[k + 1L] > 0
  ratioData <- data.frame(
    k = k[keep],
    freq_k = as.numeric(counts[k[keep]]),
    freq_k1 = as.numeric(counts[k[keep] + 1L]),
    stringsAsFactors = FALSE
  )

  ratioData$ratio <- ratioData$freq_k1 / ratioData$freq_k
  ratioData$logRatio <- log(ratioData$ratio)
  ratioData$weight <- 1 / (1 / ratioData$freq_k + 1 / ratioData$freq_k1)
  ratioData$oneInflation <- as.numeric(ratioData$k == 1)
  ratioData
}

singleRratioRegInternalModelFormula <- function(ratioFormula,
                                                response = "logRatio",
                                                addOneInflation = FALSE) {
  tt <- stats::terms(ratioFormula)
  labels <- attr(tt, "term.labels")
  if (isTRUE(addOneInflation)) {
    labels <- c(labels, "oneInflation")
  }
  stats::reformulate(
    termlabels = labels,
    response = response,
    intercept = attr(tt, "intercept")
  )
}

singleRratioRegInternalFitModel <- function(ratioData, ratioFormula, modelName) {
  baseFormula <- singleRratioRegInternalModelFormula(ratioFormula)
  fitFormula <- singleRratioRegInternalModelFormula(
    ratioFormula = ratioFormula,
    addOneInflation = identical(modelName, "M1")
  )

  if (identical(modelName, "M1") && !any(ratioData$oneInflation == 1)) {
    return(list(
      model = modelName,
      available = FALSE,
      message = "Model M1 requires an observed ratio for k = 1."
    ))
  }

  design <- stats::model.matrix(fitFormula, data = ratioData)
  if (nrow(ratioData) < ncol(design)) {
    return(list(
      model = modelName,
      available = FALSE,
      message = "Not enough usable ratios to identify the requested ratio-regression model."
    ))
  }

  if (qr(design)$rank < ncol(design)) {
    return(list(
      model = modelName,
      available = FALSE,
      message = "The ratio regression design matrix is rank deficient."
    ))
  }

  fit <- stats::lm(
    formula = fitFormula,
    data = ratioData,
    weights = weight
  )

  coefficients <- stats::coef(fit)
  if (anyNA(coefficients)) {
    return(list(
      model = modelName,
      available = FALSE,
      message = "The requested model produced aliased coefficients."
    ))
  }

  nRatio <- nrow(ratioData)
  sse <- sum(ratioData$weight * stats::residuals(fit) ^ 2)
  sigma2 <- max(sse / nRatio, .Machine$double.eps)
  logL <- -0.5 * nRatio * (log(2 * pi) + 1 + log(sigma2)) +
    0.5 * sum(log(ratioData$weight))
  nPar <- fit$rank + 1L

  list(
    model = modelName,
    available = TRUE,
    fit = fit,
    formula = fitFormula,
    baseFormula = baseFormula,
    coefficients = coefficients,
    baseCoefficients = if (identical(modelName, "M1")) {
      coefficients[names(coefficients) != "oneInflation"]
    } else {
      coefficients
    },
    rank = fit$rank,
    sigma2 = sigma2,
    logL = logL,
    AIC = -2 * logL + 2 * nPar,
    BIC = -2 * logL + log(nRatio) * nPar,
    message = NULL
  )
}

singleRratioRegInternalPredictLinear <- function(fitInfo,
                                                 newdata,
                                                 component = c("full", "base")) {
  component <- match.arg(component)
  if (!isTRUE(fitInfo$available)) {
    stop("Cannot predict from an unavailable ratio-regression fit.")
  }

  if (identical(component, "full")) {
    stats::predict(fitInfo$fit, newdata = newdata)
  } else {
    X <- stats::model.matrix(
      stats::delete.response(stats::terms(fitInfo$baseFormula)),
      data = newdata
    )
    drop(X %*% fitInfo$baseCoefficients)
  }
}

singleRratioRegInternalLogSumExp <- function(x) {
  xmax <- max(x)
  xmax + log(sum(exp(x - xmax)))
}

singleRratioRegInternalBaseProbabilities <- function(fitInfo, supportMax) {
  ratioGrid <- data.frame(
    k = 0:(supportMax - 1L),
    oneInflation = as.numeric(0:(supportMax - 1L) == 1L)
  )
  logRatios <- singleRratioRegInternalPredictLinear(
    fitInfo = fitInfo,
    newdata = ratioGrid,
    component = "base"
  )
  logTerms <- c(0, cumsum(logRatios))
  logDenom <- singleRratioRegInternalLogSumExp(logTerms)
  probabilities <- exp(logTerms - logDenom)
  names(probabilities) <- as.character(0:supportMax)
  list(
    probabilities = probabilities,
    baseRatios = exp(logRatios)
  )
}

singleRratioRegInternalPlotX <- function(baseFormula, newdata) {
  X <- stats::model.matrix(
    stats::delete.response(stats::terms(baseFormula)),
    data = newdata
  )
  keep <- colnames(X) != "(Intercept)"
  if (sum(keep) == 1L) {
    x <- as.numeric(X[, keep])
    attr(x, "label") <- colnames(X)[keep]
    return(x)
  }

  x <- newdata$k
  attr(x, "label") <- "k"
  x
}

singleRratioRegInternalPopulationEstimates <- function(fitInfo, counts) {
  baseInfo <- singleRratioRegInternalBaseProbabilities(
    fitInfo = fitInfo,
    supportMax = length(counts)
  )
  probs <- baseInfo$probabilities
  p0 <- probs[1]
  p1 <- if (length(probs) >= 2L) probs[2] else 0
  f1 <- counts[1]
  f2 <- if (length(counts) >= 2L) counts[2] else 0
  sizeObserved <- sum(counts)

  hiddenHT <- if (identical(fitInfo$model, "M1")) {
    (sizeObserved - f1) * p0 / max(1 - p0 - p1, .Machine$double.eps)
  } else {
    sizeObserved * p0 / max(1 - p0, .Machine$double.eps)
  }

  hiddenSM <- f2 / max(prod(baseInfo$baseRatios[1:2]), .Machine$double.eps)

  list(
    HT = sizeObserved + hiddenHT,
    SM = sizeObserved + hiddenSM,
    baseProbabilities = probs,
    truncatedProbabilities = probs[-1] / max(1 - p0, .Machine$double.eps),
    baseRatios = baseInfo$baseRatios
  )
}

singleRratioRegInternalMakePopResult <- function(pointEstimate,
                                                 boot = NULL,
                                                 alpha = 0.05) {
  if (is.numeric(boot) && length(boot) > 1L) {
    quant <- stats::quantile(
      x = boot,
      probs = c(alpha / 2, 1 - alpha / 2),
      names = FALSE
    )
    variance <- stats::var(boot)
    confidenceInterval <- data.frame(
      lowerBound = quant[1],
      upperBound = quant[2]
    )
  } else {
    variance <- NA_real_
    confidenceInterval <- data.frame(
      lowerBound = NA_real_,
      upperBound = NA_real_
    )
    boot <- NULL
  }

  structure(
    list(
      pointEstimate = unname(pointEstimate),
      variance = variance,
      confidenceInterval = confidenceInterval,
      boot = boot,
      control = list(alpha = alpha)
    ),
    class = "popSizeEstResults"
  )
}

singleRratioRegInternalDrawSize <- function(x) {
  base <- floor(x)
  base + stats::rbinom(1, 1, x - base)
}

singleRratioRegInternalBootstrapSample <- function(object) {
  if (identical(object$requestedModel, "auto")) {
    probs <- object$counts / sum(object$counts)
    size <- singleRratioRegInternalDrawSize(sum(object$counts))
    sampled <- as.numeric(stats::rmultinom(1, size = size, prob = probs))
    names(sampled) <- names(object$counts)
    return(sampled)
  }

  pointEstimate <- object$populationSizeAll[[object$primaryEstimator]]$pointEstimate
  f0hat <- max(pointEstimate - sum(object$counts), 0)
  probs <- c(f0hat, object$counts)
  probs <- probs / sum(probs)
  size <- singleRratioRegInternalDrawSize(pointEstimate)
  sampled <- as.numeric(stats::rmultinom(1, size = size, prob = probs))
  counts <- sampled[-1]
  names(counts) <- names(object$counts)
  counts
}

singleRratioRegInternalBootstrap <- function(object, B, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  bootValues <- matrix(
    NA_real_,
    nrow = 0,
    ncol = 2,
    dimnames = list(NULL, c("HT", "SM"))
  )
  attempts <- 0L
  maxAttempts <- max(B * 10L, B + 25L)

  while (nrow(bootValues) < B && attempts < maxAttempts) {
    attempts <- attempts + 1L
    sampledCounts <- singleRratioRegInternalBootstrapSample(object)
    bootFit <- tryCatch(
      singleRratioRegInternalFitCounts(
        counts = sampledCounts,
        formula = object$formula,
        ratioFormula = object$ratioFormula,
        call = object$call,
        requestedModel = object$requestedModel,
        criterion = object$criterion,
        primaryEstimator = object$primaryEstimator,
        confintType = "none",
        B = 0L,
        seed = NULL
      ),
      error = function(e) e
    )

    if (inherits(bootFit, "error")) {
      next
    }

    values <- c(
      HT = bootFit$populationSizeAll$HT$pointEstimate,
      SM = bootFit$populationSizeAll$SM$pointEstimate
    )

    if (all(is.finite(values))) {
      bootValues <- rbind(bootValues, values)
    }
  }

  if (nrow(bootValues) < B) {
    warning(
      "Only ", nrow(bootValues), " successful bootstrap replicates were obtained ",
      "out of the requested ", B, "."
    )
  }

  bootValues
}

singleRratioRegInternalCriteriaTable <- function(fits, selectedModel) {
  data.frame(
    model = c("M0", "M1"),
    available = c(fits$M0$available, fits$M1$available),
    logLik = c(
      if (fits$M0$available) fits$M0$logL else NA_real_,
      if (fits$M1$available) fits$M1$logL else NA_real_
    ),
    AIC = c(
      if (fits$M0$available) fits$M0$AIC else NA_real_,
      if (fits$M1$available) fits$M1$AIC else NA_real_
    ),
    BIC = c(
      if (fits$M0$available) fits$M0$BIC else NA_real_,
      if (fits$M1$available) fits$M1$BIC else NA_real_
    ),
    selected = c(selectedModel == "M0", selectedModel == "M1"),
    row.names = c("M0", "M1")
  )
}

singleRratioRegInternalFitCounts <- function(counts,
                                             formula,
                                             ratioFormula,
                                             call,
                                             requestedModel,
                                             criterion,
                                             primaryEstimator,
                                             confintType,
                                             B,
                                             seed = NULL) {
  ratioData <- singleRratioRegInternalBuildRatioData(counts)
  baseFormula <- singleRratioRegInternalModelFormula(ratioFormula)
  ratioData$xPlot <- singleRratioRegInternalPlotX(baseFormula, ratioData)

  fits <- list(
    M0 = singleRratioRegInternalFitModel(ratioData, ratioFormula, "M0"),
    M1 = singleRratioRegInternalFitModel(ratioData, ratioFormula, "M1")
  )

  selectedModel <- switch(
    requestedModel,
    M0 = {
      if (!fits$M0$available) {
        stop(fits$M0$message)
      }
      "M0"
    },
    M1 = {
      if (!fits$M1$available) {
        stop(fits$M1$message)
      }
      "M1"
    },
    auto = {
      available <- c(M0 = fits$M0$available, M1 = fits$M1$available)
      if (!any(available)) {
        stop("Neither M0 nor M1 could be fitted to the available ratios.")
      }
      if (sum(available) == 1L) {
        names(which(available))[1]
      } else {
        stats <- c(M0 = fits$M0[[criterion]], M1 = fits$M1[[criterion]])
        names(which.min(stats))[1]
      }
    }
  )

  selectedFit <- fits[[selectedModel]]
  popEstimates <- singleRratioRegInternalPopulationEstimates(selectedFit, counts)
  truncatedProbabilities <- popEstimates$truncatedProbabilities
  fittedFrequencies <- sum(counts) * truncatedProbabilities
  names(fittedFrequencies) <- names(counts)

  plotGrid <- data.frame(
    k = seq_len(length(counts) - 1L),
    oneInflation = as.numeric(seq_len(length(counts) - 1L) == 1L)
  )
  plotX <- singleRratioRegInternalPlotX(baseFormula, plotGrid)
  plotXLabel <- attr(plotX, "label")

  for (k in names(fits)) {
    if (!fits[[k]]$available) {
      next
    }

    fits[[k]]$plotGrid <- plotGrid
    fits[[k]]$plotX <- plotX
    fits[[k]]$plotPredictions <- list(
      full = singleRratioRegInternalPredictLinear(
        fitInfo = fits[[k]],
        newdata = plotGrid,
        component = "full"
      ),
      base = singleRratioRegInternalPredictLinear(
        fitInfo = fits[[k]],
        newdata = plotGrid,
        component = "base"
      )
    )
  }

  object <- structure(
    list(
      call = call,
      formula = formula,
      ratioFormula = ratioFormula,
      requestedModel = requestedModel,
      selectedModel = selectedModel,
      criterion = criterion,
      primaryEstimator = primaryEstimator,
      counts = counts,
      sizeObserved = sum(counts),
      supportMax = length(counts),
      ratioData = ratioData,
      fits = fits,
      criteria = singleRratioRegInternalCriteriaTable(fits, selectedModel),
      fittedProbabilities = popEstimates$baseProbabilities,
      fittedFrequencies = fittedFrequencies,
      plotInfo = list(
        xLabel = plotXLabel,
        xObserved = ratioData$xPlot,
        yObserved = ratioData$logRatio
      ),
      populationSize = NULL,
      populationSizeAll = NULL,
      bootstrap = NULL
    ),
    class = "singleRRatioReg"
  )

  bootstrap <- NULL
  if (identical(confintType, "bootstrap")) {
    bootstrap <- singleRratioRegInternalBootstrap(
      object = object,
      B = as.integer(B),
      seed = seed
    )
  }

  object$bootstrap <- bootstrap
  object$populationSizeAll <- list(
    HT = singleRratioRegInternalMakePopResult(
      pointEstimate = popEstimates$HT,
      boot = if (is.numeric(bootstrap)) bootstrap[, "HT"] else NULL
    ),
    SM = singleRratioRegInternalMakePopResult(
      pointEstimate = popEstimates$SM,
      boot = if (is.numeric(bootstrap)) bootstrap[, "SM"] else NULL
    )
  )
  object$populationSize <- object$populationSizeAll[[primaryEstimator]]
  object
}

#' @method coef singleRRatioReg
#' @exportS3Method
coef.singleRRatioReg <- function(object,
                                 model = c("selected", "M0", "M1"),
                                 component = c("full", "base"),
                                 ...) {
  model <- match.arg(model)
  component <- match.arg(component)
  if (identical(model, "selected")) {
    model <- object$selectedModel
  }

  fitInfo <- object$fits[[model]]
  if (!isTRUE(fitInfo$available)) {
    stop("Requested model is not available in this ratioReg() fit.")
  }

  switch(
    component,
    full = fitInfo$coefficients,
    base = fitInfo$baseCoefficients
  )
}

#' @method logLik singleRRatioReg
#' @importFrom stats logLik
#' @exportS3Method
logLik.singleRRatioReg <- function(object,
                                   model = c("selected", "M0", "M1"),
                                   ...) {
  model <- match.arg(model)
  if (identical(model, "selected")) {
    model <- object$selectedModel
  }

  fitInfo <- object$fits[[model]]
  if (!isTRUE(fitInfo$available)) {
    stop("Requested model is not available in this ratioReg() fit.")
  }

  val <- fitInfo$logL
  attr(val, "nobs") <- nrow(object$ratioData)
  attr(val, "df") <- fitInfo$rank + 1L
  class(val) <- "logLik"
  val
}

#' @method extractAIC singleRRatioReg
#' @importFrom stats extractAIC
#' @exportS3Method
extractAIC.singleRRatioReg <- function(fit,
                                       scale,
                                       k = 2,
                                       model = c("selected", "M0", "M1"),
                                       ...) {
  model <- match.arg(model)
  if (identical(model, "selected")) {
    model <- fit$selectedModel
  }

  fitInfo <- fit$fits[[model]]
  if (!isTRUE(fitInfo$available)) {
    stop("Requested model is not available in this ratioReg() fit.")
  }

  -2 * fitInfo$logL + k * (fitInfo$rank + 1L)
}

#' @method nobs singleRRatioReg
#' @importFrom stats nobs
#' @exportS3Method
nobs.singleRRatioReg <- function(object, ...) {
  nrow(object$ratioData)
}

#' @method popSizeEst singleRRatioReg
#' @rdname popSizeEst
#' @exportS3Method
popSizeEst.singleRRatioReg <- function(object,
                                       estimator = c("primary", "HT", "SM"),
                                       ...) {
  estimator <- match.arg(estimator)
  if (identical(estimator, "primary")) {
    estimator <- object$primaryEstimator
  }
  object$populationSizeAll[[estimator]]
}

#' @method confint singleRRatioReg
#' @importFrom stats confint
#' @exportS3Method
confint.singleRRatioReg <- function(object,
                                    parm,
                                    level = 0.95,
                                    estimator = c("primary", "HT", "SM"),
                                    ...) {
  if (!missing(level) && (length(level) != 1L || !is.finite(level))) {
    stop("Argument level must be a single finite number.")
  }

  estimator <- match.arg(estimator)
  if (identical(estimator, "primary")) {
    estimator <- object$primaryEstimator
  }

  boot <- object$populationSizeAll[[estimator]]$boot
  if (is.numeric(boot) && length(boot) > 1L) {
    alpha <- 1 - level
    return(
      matrix(
        stats::quantile(
          x = boot,
          probs = c(alpha / 2, 1 - alpha / 2),
          names = FALSE
        ),
        nrow = 1,
        dimnames = list(NULL, c("lowerBound", "upperBound"))
      )
    )
  }

  as.matrix(object$populationSizeAll[[estimator]]$confidenceInterval)
}

#' @method print singleRRatioReg
#' @exportS3Method
print.singleRRatioReg <- function(x, ...) {
  cat("Call: ")
  print(x$call)
  cat(
    "\nSelected ratio regression model:", x$selectedModel,
    "\nSelection criterion:", x$criterion,
    "\nPrimary population estimator:", x$primaryEstimator,
    "\nObserved population size:", signif(x$sizeObserved),
    "\n"
  )
  cat("\nCoefficients:\n")
  print(stats::coef(x))
  cat("\nPopulation size estimation results:\n")
  print(popSizeEst(x))

  alternative <- setdiff(c("HT", "SM"), x$primaryEstimator)
  cat("\nAlternative estimator point estimate:", alternative, "=",
      signif(popSizeEst(x, estimator = alternative)$pointEstimate), "\n")

  invisible(x)
}

#' @title Summary Statistics for Ratio Regression Fits
#'
#' @description
#' Summarizes a fitted `singleRRatioReg` object, including model-selection
#' criteria and population-size estimates.
#'
#' @param object object of class `singleRRatioReg`.
#' @param ... currently not used.
#'
#' @return An object of class `summarysingleRRatioReg`.
#' @method summary singleRRatioReg
#' @exportS3Method
summary.singleRRatioReg <- function(object, ...) {
  estimates <- data.frame(
    estimator = c("HT", "SM"),
    pointEstimate = c(
      object$populationSizeAll$HT$pointEstimate,
      object$populationSizeAll$SM$pointEstimate
    ),
    variance = c(
      object$populationSizeAll$HT$variance,
      object$populationSizeAll$SM$variance
    ),
    lowerBound = c(
      object$populationSizeAll$HT$confidenceInterval$lowerBound,
      object$populationSizeAll$SM$confidenceInterval$lowerBound
    ),
    upperBound = c(
      object$populationSizeAll$HT$confidenceInterval$upperBound,
      object$populationSizeAll$SM$confidenceInterval$upperBound
    )
  )

  structure(
    list(
      call = object$call,
      selectedModel = object$selectedModel,
      criterion = object$criterion,
      primaryEstimator = object$primaryEstimator,
      sizeObserved = object$sizeObserved,
      coefficients = lapply(
        object$fits[c("M0", "M1")],
        function(x) if (isTRUE(x$available)) x$coefficients else NULL
      ),
      criteria = object$criteria,
      estimates = estimates
    ),
    class = "summarysingleRRatioReg"
  )
}

#' @method print summarysingleRRatioReg
#' @exportS3Method
print.summarysingleRRatioReg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n")
  print(x$call)
  cat(
    "\nSelected ratio regression model:", x$selectedModel,
    "\nSelection criterion:", x$criterion,
    "\nPrimary population estimator:", x$primaryEstimator,
    "\nObserved population size:", signif(x$sizeObserved, digits),
    "\n\nModel selection table:\n",
    sep = " "
  )
  print(x$criteria, digits = digits)

  for (nm in names(x$coefficients)) {
    cat("\nCoefficients for", nm, ":\n")
    if (is.null(x$coefficients[[nm]])) {
      cat("Model not available.\n")
    } else {
      print(round(x$coefficients[[nm]], digits = digits))
    }
  }

  cat("\nPopulation size estimates:\n")
  toPrint <- x$estimates
  isNum <- vapply(toPrint, is.numeric, logical(1))
  toPrint[isNum] <- lapply(toPrint[isNum], round, digits = digits)
  print(toPrint, row.names = FALSE)
  invisible(x)
}

#' @title Diagnostic Plots for Ratio Regression Fits
#'
#' @description
#' Plot methods for `singleRRatioReg` objects.
#'
#' @param x object of class `singleRRatioReg`.
#' @param plotType one of `"ratio"`, `"bootHist"`, or `"freqFit"`.
#' @param ... additional arguments passed to the underlying graphics calls.
#'
#' @return No return value, called for side effects.
#' @method plot singleRRatioReg
#' @export
plot.singleRRatioReg <- function(x,
                                 plotType = c("ratio", "bootHist", "freqFit"),
                                 ...) {
  plotType <- match.arg(plotType)

  switch(
    plotType,
    ratio = {
      graphics::plot(
        x = x$plotInfo$xObserved,
        y = x$plotInfo$yObserved,
        xlab = x$plotInfo$xLabel,
        ylab = "log(f[k+1] / f[k])",
        main = "Ratio regression plot",
        pch = 19,
        ...
      )

      if (isTRUE(x$fits$M0$available)) {
        graphics::lines(
          x = x$fits$M0$plotX,
          y = x$fits$M0$plotPredictions$full,
          col = 2,
          lwd = if (x$selectedModel == "M0") 2.5 else 1.5,
          lty = 2
        )
      }

      if (isTRUE(x$fits$M1$available)) {
        graphics::lines(
          x = x$fits$M1$plotX,
          y = x$fits$M1$plotPredictions$full,
          col = 4,
          lwd = if (x$selectedModel == "M1") 2.5 else 1.5,
          lty = 1
        )
      }

      graphics::legend(
        "topright",
        legend = c("Observed", "M0", "M1"),
        col = c(1, 2, 4),
        pch = c(19, NA, NA),
        lty = c(NA, 2, 1),
        bty = "n"
      )
    },
    bootHist = {
      boot <- popSizeEst(x)$boot
      if (!is.numeric(boot)) {
        stop("No bootstrap results were stored for the primary estimator.")
      }

      graphics::hist(
        boot,
        main = "Bootstrap population size estimates",
        xlab = expression(hat(N)),
        ...
      )
    },
    freqFit = {
      xx <- seq_len(length(x$counts))
      graphics::matplot(
        x = xx,
        y = cbind(observed = x$counts, fitted = x$fittedFrequencies),
        type = "o",
        pch = c(19, 17),
        lty = c(1, 2),
        col = c(1, 2),
        xlab = "Counts",
        ylab = "Frequency",
        main = "Observed and fitted frequencies",
        ...
      )
      graphics::legend(
        "topright",
        legend = c("Observed", "Fitted"),
        col = c(1, 2),
        pch = c(19, 17),
        lty = c(1, 2),
        bty = "n"
      )
    }
  )

  invisible(x)
}
