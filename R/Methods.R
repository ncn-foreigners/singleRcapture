#' @title Summary statistics for model of singleRStaticCountData class.
#' 
#' @description A \code{summary} method for \code{singleRStaticCountData} class, works 
#' analogically to \code{summary.glm} but includes population size estimation 
#' results. If any additional statistics, such as confidence intervals for
#' coefficients or coefficient correlation, are specified they will be printed.
#' 
#' @param object object of singleRStaticCountData class.
#' @param test type of test for significance of parameters \code{"t"} for t-test 
#' and \code{"z"} for normal approximation of students t distribution, by 
#' default \code{"z"} is used if there are more than 30 degrees of freedom
#' and \code{"t"} is used in other cases.
#' @param resType type of residuals to summarize any value that is allowed in
#' \code{residuals.singleRStaticCountData} except for \code{"all"} is allowed. By default 
#' pearson residuals are used.
#' @param correlation logical value indicating whether correlation matrix should
#' be computed from covariance matrix by default \code{FALSE}.
#' @param confint logical value indicating whether confidence intervals for
#' regression parameters should be constructed. By default \code{FALSE}.
#' @param cov covariance matrix corresponding to regression parameters. 
#' It is possible to give \code{cov} argument as a function of \code{object}.
#' If not specified it will be constructed using \code{vcov.singleRStaticCountData} method.
#' (i.e using Cramer-Rao lower bound)
#' @param popSizeEst a \code{popSizeEstResults} class object.
#' If not specified population size estimation results will be drawn from
#' \code{object}. If any post-hoc procedures, such as sandwich covariance matrix 
#' estimation or bias reduction, were taken it is possible to include them in 
#' population size estimation results by calling \code{redoPopEstimation}.
#' @param ... additional optional arguments passed to the following functions:
#' \itemize{
#' \item \code{vcov.singleRStaticCountData} -- if no \code{cov} argument was provided.
#' \item \code{cov} -- if \code{cov} parameter specified at call was a function.
#' \item \code{confint.singleRStaticCountData} -- if \code{confint} parameter was set to \code{TRUE} at function call.
#' In particular it is possible to set confidence level in \code{...}.
#' }
#' 
#' @return An object of \code{summarysingleRStaticCountData} class containing:
#' \itemize{
#' \item \code{call} -- A call which created \code{object}.
#' \item \code{coefficients} -- A dataframe with estimated regression coefficients
#' and their summary statistics such as standard error Wald test statistic and
#' p value for Wald test.
#' \item \code{residuals} -- A vector of residuals of type specified at call.
#' \item \code{aic} -- Akaike's information criterion.
#' \item \code{bic} -- Bayesian (Schwarz's) information criterion.
#' \item \code{iter} -- Number of iterations taken in fitting regression.
#' \item \code{logL} -- Logarithm of likelihood function evaluated at coefficients.
#' \item \code{deviance} -- Residual deviance.
#' \item \code{populationSize} -- Object with population size estimation results.
#' \item \code{dfResidual} -- Residual degrees of freedom.
#' \item \code{sizeObserved} -- Size of observed population.
#' \item \code{correlation} -- Correlation matrix if \code{correlation} parameter was set to \code{TRUE}
#' \item \code{test} -- Type of statistical test performed.
#' \item \code{model} -- Family class object specified in call for \code{object}.
#' \item \code{skew} -- If bootstrap sample was saved contains estimate of skewness.
#' }
#' @method summary singleRStaticCountData
#' @importFrom stats pt
#' @importFrom stats coef
#' @importFrom stats sd
#' @seealso [redoPopEstimation()] [stats::summary.glm()]
#' @exportS3Method 
summary.singleRStaticCountData <- function(object, 
                                           test = c("t", "z"), 
                                           resType = "pearson", 
                                           correlation = FALSE, 
                                           confint = FALSE, 
                                           cov, 
                                           popSizeEst, 
                                           ...) {
  if (resType == "all") {stop("Can't use 'resType = all' in summary.singleRStaticCountData method, if you wish to obtain all aviable types of residuals call residuals.singleRStaticCountData method directly.")}
  dfResidual <- object$dfResidual
  if (missing(test)) {if (dfResidual > 30) test <- "z" else test <- "t"}
  if (missing(cov)) {
    cov <- vcov(object, ...)
  } else if (is.function(cov)) {
    cov <- cov(object, ...)
  }
  
  pers <- residuals(object, type = resType)
  
  cf <- stats::coef(object)
  se <- sqrt(diag(cov))
  
  wValues <- cf / se
  
  pValues <- switch (test,
    "t" = 2 *    stats::pt(q = -abs(wValues), df = dfResidual),
    "z" = 2 * stats::pnorm(q = abs(wValues), lower.tail = FALSE)
  )
  
  crr <- if (isFALSE(correlation)) {NULL} else {cov / outer(se, se)}
  
  if(isTRUE(correlation)) {rownames(crr) <- colnames(crr) <- names(cf)}
  
  cnfint <- if(isTRUE(confint)) {confint(object, ...)} else {NULL}
  
  if (is.numeric(object$populationSize$boot)) {
    n <- length(object$populationSize$boot)
    m <- sum((object$populationSize$boot - mean(object$populationSize$boot)) ^ 3) / n
    s <- sd(object$populationSize$boot)
    skew <- m / (s ^ 3)
  } else {
    skew <- NULL
  }
  
  ob <- data.frame(cf, se, wValues, pValues)
  
  colnames(ob) <- switch(test,
    "t" = c("Estimate", "Std. Error", "t value", "P(>|t|)"),
    "z" = c("Estimate", "Std. Error", "z value", "P(>|z|)")
  )
  if (isTRUE(confint)) {
    ob[, 4] <- cnfint[, 1]
    ob[, 5] <- cnfint[, 2]
    ob[, 6] <- pValues
    colnames(ob)[4:6] <- c(colnames(cnfint), colnames(ob)[4])
  }
  
  cf <- ob
  structure(
    list(
      call = object$call,
      coefficients = cf,
      residuals = pers,
      aic = AIC(object, ...),
      bic = BIC(object, ...),
      iter = object$iter,
      logL = object$logL,
      deviance = object$deviance,
      populationSize = if (missing(popSizeEst)) object$populationSize else popSizeEst,
      dfResidual = dfResidual,
      sizeObserved = object$sizeObserved,
      correlation = crr,
      test = test,
      model = object$model,
      skew = skew
    ),
    class = "summarysingleRStaticCountData"
  )
}

#' Predict method for \code{singleRStaticCountData} class
#'
#' \loadmathjax
#' 
#' @description
#' A method for \code{predict} function, works analogous to \code{predict.glm}
#' but gives the possibility to get standard errors of 
#' mean/distribution parameters and directly get pop size estimates for new data.
#' 
#' 
#' @param object an object of \code{singleRStaticCountData} class.
#' @param newdata an optional \code{data.frame} containing new data.
#' @param type the type of prediction required, possible values are:
#' \itemize{
#'   \item \code{"response"}-- For matrix containing estimated distributions
#'   parameters.
#'   \item \code{"link"}    -- For matrix of linear predictors.
#'   \item \code{"mean"}    -- For fitted values of both \mjseqn{Y} and
#'   \mjseqn{Y|Y>0}.
#'   \item \code{"contr"}   -- For inverse probability weights (here named for 
#'   observation contribution to population size estimate).
#'   \item \code{"popSize"} -- For population size estimation. Note
#'   this results in a call to \code{redoPopEstimation} and it is
#'   usually better to call this function directly.
#' } by default set to \code{"response"}.
#' @param se.fit a logical value indicating whether standard errors should be 
#' computed. Only matters for \code{type} in \code{"response", "mean", "link"}.
#' @param na.action does nothing yet.
#' @param weights optional vector of weights for \code{type} in \code{"contr", "popSize"}.
#' @param cov optional matrix or function or character specifying either
#' a covariance matrix or a function to compute that covariance matrix.
#' By default \code{vcov.singleRStaticCountData} can be set to e.g. \code{vcovHC}.
#' @param ... arguments passed to other functions, for now this only affects
#' \code{vcov.singleRStaticCountData} method and \code{cov} function.
#' 
#' @details Standard errors are computed with assumption of regression
#' coefficients being asymptotically normally distributed, if this assumption
#' holds then each of linear predictors i.e. each row of
#' \mjseqn{\boldsymbol{\eta}=\boldsymbol{X}_{vlm}\boldsymbol{\beta}}
#' is asymptotically normally distributed and their variances are expressed by
#' well known formula. The mean \mjseqn{\mu} and distribution parameters
#' are then differentiable functions of asymptotically normally distributed
#' variables and therefore their variances can be computed using (multivariate)
#' delta method.
#'
#' @return Depending on \code{type} argument if one of \code{"response", "link", "mean"}
#' a matrix with fitted values and possibly standard errors if \code{se.fit}
#' argument was set to \code{TRUE}, if \code{type} was set to \code{"contr"}
#' a vector with inverses of probabilities, finally for \code{"popSize"}
#' an object of class \code{popSizeEstResults} with its own methods containing
#' population size estimation results.
#' 
#' @method predict singleRStaticCountData
#' @seealso [redoPopEstimation()] [stats::summary.glm()] [estimatePopsize()]
#' @exportS3Method
predict.singleRStaticCountData <- function(object,
                                           newdata,
                                           type = c("response", "link", 
                                                    "mean", "popSize", 
                                                    "contr"),
                                           se.fit = FALSE,
                                           na.action = NULL,
                                           weights,
                                           cov,
                                           ...) {
  type <- match.arg(type)
  if (missing(weights)) {
    if (missing(newdata)) {
      weights <- object$priorWeights
    } else {
      weights <- rep(1, NROW(newdata))
    }
  }
  
  if (missing(cov)) {
    cov <- vcov
  }
  
  if (is.character(cov)) {
    cov <- get(cov, mode = "function", envir = parent.frame())
  }
  if (is.function(cov)) {
    cov <- cov(object, ...)
  }
  

  if (missing(newdata)) {
    Xvlm <- model.matrix(
      object, type = "vlm"
    )
    
    eta <- object$linearPredictors
    
    res <- switch (type,
      response = as.data.frame(
        lapply(1:length(family(object)$etaNames), FUN = function(x) {
          family(object)$links[[x]](
            eta[, x],
            inverse = TRUE
          )
      }), col.names = family(object)$etaNames),
      link     = eta,
      mean     = fitted(object, "all"),
      popSize  = popSizeEst(object, ...),
      contr    = family(object)$pointEst(
        pw     = weights,
        eta    = eta,
        contr  = TRUE,
        y      = model.response(model.frame(object))
      )
    )
  } else {
    mf <- model.frame(
      object, data = newdata
    )
    
    Xvlm <- model.matrix(
      object, type = "vlm", data = newdata
    )
    
    eta <- matrix(
      as.matrix(Xvlm) %*% stats::coef(object), 
      ncol = length(family(object)$etaNames),
      dimnames = list(
        rownames(Xvlm),
        family(object)$etaNames
      )
    )

    res <- switch (type,
      response = as.data.frame(
        lapply(1:length(family(object)$etaNames), FUN = function(x) {
        family(object)$links[[x]](
          eta[, x],
          inverse = TRUE
        )
      }), col.names = family(object)$etaNames),
      link = eta,
      mean = data.frame(
        "truncated"    = family(object)$mu.eta(eta = eta),
        "nontruncated" = family(object)$mu.eta(eta = eta, type = "nontrunc")
      ),
      popSize = redoPopEstimation(
        object = object, newdata = newdata,
        weights = if (missing(weights)) rep(1, NROW(mf)) else weights,
        cov = cov,
        ...
      ),
      contr = family(object)$pointEst(
        pw    = weights,
        eta   = eta,
        contr = TRUE,
        y     = model.response(mf) 
      )
    )
  }
  
  if (isTRUE(se.fit)) {
    cov <- vcov(object, ...)
    
    #### beta is asymptotically normal and each eta is just a linear
    #### combination of beta so it is also asymptotically normal
    #### and its variance is computed by quadratic form
    if (type != "mean") {
      se <- matrix(
        sapply(1:NROW(Xvlm), function(x) {
          (t(Xvlm[x,]) %*% cov %*% Xvlm[x,]) ^ .5
        }),
        ncol = length(family(object)$etaNames),
        dimnames = list(
          rownames(eta),
          paste0("se:", family(object)$etaNames)
        )
      )
      
      if (type == "link") 
        res <- cbind(res, se)
      
      ## since beta is asymptotically normal we use delta metod for
      ## geting standard errors of invlink(eta)
      
      if (type == "response") {
        se <- matrix(
          sapply(
            1:length(family(object)$etaNames),
            function (x) {
              se[, x, drop = FALSE] *
                abs(family(object)$links[[x]](
                  eta[, x],
                  inverse = TRUE,
                  deriv = 1
                ))
            }
          ),
          ncol = length(family(object)$etaNames),
          dimnames = dimnames(se)
        )
        res <- cbind(res,se)
      }
    } else if (type == "mean") {
      ## covariance matrix foe each row of linear predictors
      auxVec <- c(0, cumsum(rep(
        NROW(eta), 
        length.out = length(family(object)$etaNames) - 1
      )))
      
      se <- lapply(
        1:NROW(eta),
        function (x) {
          Xvlm[x + auxVec, , drop = FALSE] %*% cov %*%
            t(Xvlm[x + auxVec, , drop = FALSE])
        }
      )
      
      derivMu <- list(
        family(object)$mu.eta(eta, type = "nontrunc", deriv = 1),
        family(object)$mu.eta(eta, type = "trunc", deriv = 1)
      )
      
      res <- data.frame(
        res,
        "se:truncated" = sapply(
          1:NROW(eta),
          function (x) {
            (derivMu[[2]][x, , drop = FALSE] %*% se[[x]] %*%
               t(derivMu[[2]][x, , drop = FALSE])) ^ .5
          }
        ),
        "se:nontruncated" = sapply(
          1:NROW(eta), 
          function (x) {
            (derivMu[[1]][x, , drop = FALSE] %*% se[[x]] %*%
               t(derivMu[[1]][x, , drop = FALSE])) ^ .5
          }
        )
      )
    }
  }

  res
}

#' @title Updating population size estimation results.
#'
#' @description A function that applies all post-hoc procedures that were taken
#' (such as heteroscedastic consistent covariance matrix estimation or bias
#' reduction) to population size estimation and standard error estimation.
#' 
#' @param object object for which update of population size estimation results will be done.
#' @param newdata optional \code{data.frame} with new data for pop size estimation.
#' @param cov an updated covariance matrix estimate.
#' @param coef optional vector of coefficients of regression on which to base 
#' population size estimation. If missing it is set to \code{coef(object)}.
#' @param weights optional vector of weights to use in population size estimation. 
#' @param control similar to \code{controlPopVar} in [estimatePopsize()].
#' If missing set to controls provided on call to \code{object}.
#' @param popVar similar to \code{popVar} in [estimatePopsize()].
#' If missing set to \code{"analytic"}.
#' @param offset offset argument for new data
#' @param weightsAsCounts for \code{singleRStaticCountData} method used to specify
#' whether weights should be treated as number of occurences for rows in data
#' @param ... additional optional arguments, currently not used in \code{singleRStaticCountData} class method.
#'
#' @return An object of class \code{popSizeEstResults} containing updated 
#' population size estimation results.
#' 
#' @examples
#' # Create simple model
#' Model <- estimatePopsize(
#'   formula = capture ~ nation + gender, 
#'   data = netherlandsimmigrant, 
#'   model = ztpoisson, 
#'   method = "IRLS"
#' )
#' # Apply heteroscedasticity consistent covariance matrix estimation
#' require(sandwich)
#' cov <- vcovHC(Model, type = "HC3")
#' summary(Model, cov = cov,
#' popSizeEst = redoPopEstimation(Model, cov = cov))
#' # Compare to results with usual covariance matrix estimation
#' summary(Model)
#' 
#' ## get confidence interval with larger significance level
#' redoPopEstimation(Model, control = controlPopVar(alpha = .000001))
#' @export
redoPopEstimation <- function(object, newdata, ...) {
  UseMethod("redoPopEstimation")
}

#' @title Extract population size estimation results.
#' 
#' @description An extractor function with \code{singleRStaticCountData} method for extracting
#' important information regarding pop size estimate.
#'
#' @param object object with population size estimates.
#' @param ... additional optional arguments, currently not used in \code{singleRStaticCountData} class method. 
#'
#' @return An object of class \code{popSizeEstResults} containing population size estimation results.
#' @export
popSizeEst <- function(object, ...) {
  UseMethod("popSizeEst")
}

#' @title Estimate size of sub populations.
#' 
#' @description A function that estimates sizes of specific sub populations 
#' based on a capture-recapture model for the whole population.
#'
#' @param object an object on which the population size estimates should be based
#' in \code{singleRcapture} package this is a fitter \code{singleRStaticCountData} class object.
#' @param stratas a specification of sub populations either by:
#' \itemize{
#' \item formula -- a formula to be applied to \code{model.frame} extracted from
#' the object .
#' \item Logical vector with number of entries equal to number of rows in the dataset.
#' \item A (named) list where each element is a logical vector, names of the list
#' will be used to specify names variable in returned object.
#' \item Vector of names of explanatory variables. For \code{singleRStaticCountData} method
#' for this function this specification of \code{stratas} parameter will
#' result in every level of explanatory variable having its own sub population
#' for each variable specified.
#' \item If no value was provided the \code{singleRStaticCountData} method for this function 
#' will itself create sub populations based on levels of factor variables
#' in \code{model.frame}.
#' }
#' @param cov for \code{singleRStaticCountData} method an estimate of variance-covariance matrix
#' for estimate of regression parameters. It is possible to pass a function
#' such as for example \code{sandwich::vcovHC} which will be called as:
#' \code{foo(object, ...)} and a user may specify additional arguments of a 
#' function in \code{...} argument. If not provided an estimate for covariance
#' matrix will be set by calling appropriate \code{vcov} method.
#' @param alpha significance level for confidence intervals --
#' Either a single numeric value or a vector of length equal to number of 
#' sub populations specified in \code{stratas}. 
#' If missing it is set to \code{.05} in \code{singleRStaticCountData} method.
#' @param ... a vector of arguments to be passed to other functions.
#' For \code{singleRStaticCountData} method for this functions arguments in \code{...} are 
#' passed to either \code{cov} if argument provided was a function or 
#' \code{vcov} if \code{cov} argument was missing at call.
#' 
#' \loadmathjax
#' @details In single source capture-recapture models the most frequently used
#' estimate for population size is Horvitz-Thompson type estimate:
#' 
#' \mjsdeqn{\hat{N} = \sum_{k=1}^{N}\frac{I_{k}}{\mathbb{P}(Y_{k}>0)} = 
#' \sum_{k=1}^{N_{obs}}\frac{1}{1-\mathbb{P}(Y_{k}=0)}}
#'
#' where \mjseqn{I_{k}=I_{Y_{k} > 0}} are 
#' indicator variables, with value 1 if kth unit was observed at least once 
#' and 0 otherwise and the inverse probabilistic weights weights for 
#' units observed in the data \mjseqn{\tfrac{1}{\mathbb{P}(Y_{k}>0)}}
#' are estimated using fitted linear predictors.
#' 
#' The estimates for different sub populations are made by changing the
#' \mjseqn{I_{k}=I_{Y_{k} > 0}} indicator variables to 
#' refer not to the population as a whole but to the sub populations that are 
#' being considered i.e. by changing values from 1 to 0 if kth unit is not a 
#' member of sub population that is being considered at the moment.
#' 
#' The estimation of variance for these estimates and estimation of variance for
#' estimate of population size for the whole population follow the same relation
#' as the one described above.
#' 
#' @seealso [vcov.singleRStaticCountData()] [estimatePopsize()]
#'
#' @return A \code{data.frame} object with row names being the names of specified 
#' sub populations either provided or inferred.
#' @export
stratifyPopsize <- function(object, stratas, alpha, ...) {
  UseMethod("stratifyPopsize")
}

#' @title Obtain Covariance Matrix estimation.
#' 
#' @description A \code{vcov} method for \code{singleRStaticCountData} class.
#' 
#' @param object object of singleRStaticCountData class.
#' @param type type of estimate for covariance matrix for now either
#' expected (Fisher) information matrix or observed information matrix.
#' @param ... additional arguments for method functions
#' 
#' @details  Returns a estimated covariance matrix for model coefficients
#' calculated from analytic hessian or Fisher information matrix usually 
#' utilizing asymptotic effectiveness of maximum likelihood estimates.
#' Covariance type is taken from control parameter that have been provided
#' on call that created \code{object} if arguments \code{type} was not specified.
#' 
#' @method vcov singleRStaticCountData
#' @return A covariance matrix for fitted coefficients, rows and columns of which 
#' correspond to parameters returned by \code{coef} method.
#' @seealso [vcovHC.singleRStaticCountData()] [sandwich::sandwich()]
#' @exportS3Method
vcov.singleRStaticCountData <- function(object, 
                         type = c("Fisher", 
                                  "observedInform"), 
                         ...) {
  if (missing(type) ){type <- object$populationSize$control$covType}
  X <- model.matrix(object, "vlm")
  res <- switch(
    type,
    "observedInform" = solve(
      -object$model$makeMinusLogLike(y = object$y, X = X,
      weight = object$priorWeights, deriv = 2, 
      offset = object$offset)(object$coefficients),
      ...
    ),
    "Fisher" = {
      if (isTRUE(object$call$method == "IRLS")) {W <- object$weights} else {W <- object$model$Wfun(prior = object$priorWeights, eta = object$linearPredictors, y = as.numeric(model.response(model.frame(object))))};
      if (isTRUE(object$control$controlMethod$checkDiagWeights)) {
        W[, (1:length(family(object)$etaNames)) ^ 2] <- ifelse(
          W[, (1:length(family(object)$etaNames)) ^ 2] < object$control$controlMethod$weightsEpsilon, 
          object$control$controlMethod$weightsEpsilon, 
          W[, (1:length(family(object)$etaNames)) ^ 2]
        )
      };
      solve(
      singleRinternalMultiplyWeight(X = X, W = W) %*% X,
      ...
    )}
  )
  dimnames(res) <- list(names(object$coefficients), names(object$coefficients))
  res
}

#' @title Confidence Intervals for Model Parameters
#' 
#' @description A function that computes studentized confidence intervals
#' for model coefficients.
#' 
#' @param object object of singleRStaticCountData class.
#' @param parm names of parameters for which confidence intervals are to be 
#' computed, if missing all parameters will be considered.
#' @param level confidence level for intervals.
#' @param ... currently does nothing.
#' 
#' @method confint singleRStaticCountData
#' @return An object with named columns that include upper and 
#' lower limit of confidence intervals.
#' @exportS3Method
confint.singleRStaticCountData <- function(object,
                            parm, 
                            level = 0.95, 
                            ...) {
  std <- sqrt(diag(vcov.singleRStaticCountData(object)))
  if (missing(parm)) {
    coef <- object$coefficients
  } else {
    coef <- object$coefficients[parm]
    std <- std[parm]
  }
  sc <- qnorm(p = 1 - (1 - level) / 2)
  res <- data.frame(coef - sc * std, coef + sc * std)
  colnames(res) <- c(paste0(100 * (1 - level) / 2, "%"),
                     paste0(100 * (1 - (1 - level) / 2), "%"))
  res
}

# functions with documented elsewhere

#' @rdname regDiagSingleR
#' @export
dfpopsize <- function(model, ...) {
  UseMethod("dfpopsize")
}

#' @rdname regDiagSingleR
#' @method hatvalues singleRStaticCountData
#' @importFrom stats hatvalues
#' @exportS3Method 
hatvalues.singleRStaticCountData <- function(model, ...) {
  X <- model.frame.singleRStaticCountData(model, ...)
  X <- singleRinternalGetXvlmMatrix(
    X = X, 
    formulas = model$formula, 
    parNames = model$model$etaNames
  )
  if (isTRUE(model$call$method == "IRLS")) {
    W <- model$weights
  } else {
    W <- model$model$Wfun(
      prior = model$priorWeights, 
      eta   = model$linearPredictors,
      y     = if (is.null(model$y)) model.response(model.frame(model)) 
      else model$y
    )
  }
  
  mlt <- singleRinternalMultiplyWeight(X = X, W = W)
  hatvalues <- diag(X %*% solve(mlt %*% X) %*% mlt)
  hatvalues <- matrix(
    hatvalues, 
    ncol = length(model$model$etaNames), 
    dimnames = list(
      1:(length(hatvalues) / length(model$model$etaNames)), 
      model$model$etaNames
    )
  )
  hatvalues
}
#' @method dfbeta singleRStaticCountData
#' @importFrom stats dfbeta
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @rdname regDiagSingleR
#' @exportS3Method 
dfbeta.singleRStaticCountData <- function(model,
                           maxitNew = 1,
                           trace = FALSE,
                           cores = 1,
                           ...) {
  # formula method removed since it doesn't give good results will reimplement if we find better formula
  X <- model.frame.singleRStaticCountData(model, ...)
  y <- if (is.null(model$y)) stats::model.response(X) else model$y
  X <- X
  y <- y
  cf <- coef(model)
  pw <- model$priorWeights
  offset <- model$offset
  eta <- model$linearPredictors
  
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    #parallel::clusterExport(cl, c("singleRinternalGetXvlmMatrix", "cf", "y", "X", "maxitNew", "model", "pw", "offset", "eta"), envir = environment())
    
    if (isFALSE(model$control$controlModel$weightsAsCounts)) {
      res <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:NROW(X), .combine = rbind),
        ex = {
          c(cf - estimatePopsizeFit(
            control = controlMethod(
              silent = TRUE,
              maxiter = maxitNew + 1,
              ...
            ),
            y = y[-k],
            X = singleRinternalGetXvlmMatrix(
              X        = X[rownames(X) != rownames(X)[k], , drop = FALSE],
              formulas = model$formula,
              parNames = model$model$etaNames
            ),
            coefStart = cf,
            etaStart  = eta[-k, , drop = FALSE] + offset[-k, , drop = FALSE],
            family = model$model,
            priorWeights = pw[-k],
            method = if (is.null(model$call$method)) "IRLS" else model$call$method,
            offset = offset[-k, , drop = FALSE]
          )$beta)
        }
      )
    } else {
      res <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:NROW(X), .combine = rbind),
        ex = {
          if (isFALSE(pw[k] - 1 > 0)) {
            c(cf - estimatePopsizeFit(
              control = controlMethod(
                silent = TRUE,
                maxiter = maxitNew + 1,
                ...
              ),
              y = y[-k],
              X = singleRinternalGetXvlmMatrix(
                X        = X[-k, , drop = FALSE],
                formulas = model$formula,
                parNames = model$model$etaNames
              ),
              coefStart = cf,
              etaStart  = eta[-k, , drop = FALSE] + offset[-k, , drop = FALSE],
              family = model$model,
              priorWeights = pw[-k],
              method = if (is.null(model$call$method)) "IRLS" else model$call$method,
              offset = offset[-k, , drop = FALSE]
            )$beta)
          } else {
            kk <- rep(0, length(pw))
            kk[k] <- 1
            
            c(cf - estimatePopsizeFit(
              control = controlMethod(
                silent = TRUE,
                maxiter = maxitNew + 1,
                ...
              ),
              y = y,
              X = model.matrix(model, "vlm"),
              coefStart = cf,
              etaStart  = eta + offset,
              family = model$model,
              priorWeights = pw - kk,
              method = if (is.null(model$call$method)) "IRLS" else model$call$method,
              offset = offset
            )$beta)
          }
        }
      )
    }
  } else {
    res <- matrix(nrow = nrow(X), ncol = length(cf))
    
    for (k in 1:nrow(X)) {
      if (isTRUE(trace)) {
        cat("-----\nRemoving observation number: ", k, "\n", sep = "")
      }
      if (isFALSE(model$control$controlModel$weightsAsCounts)) {
        res[k, ] <- cf - estimatePopsizeFit(
          control = controlMethod(
            silent = TRUE,
            maxiter = maxitNew + 1,
            ...
          ),
          y = y[-k],
          X = singleRinternalGetXvlmMatrix(
            X        = X[rownames(X) != rownames(X)[k], , drop = FALSE],
            formulas = model$formula, 
            parNames = model$model$etaNames
          ),
          coefStart = cf,
          etaStart  = eta[-k, , drop = FALSE] + offset[-k, , drop = FALSE],
          family = model$model,
          priorWeights = pw[-k],
          method = if (is.null(model$call$method)) "IRLS" else model$call$method,
          offset = offset[-k, , drop = FALSE]
        )$beta
      } else {
        if (isFALSE(pw[k] - 1 > 0)) {
          res[k, ] <- cf - estimatePopsizeFit(
            control = controlMethod(
              silent = TRUE,
              maxiter = maxitNew + 1,
              ...
            ),
            y = y[-k],
            X = singleRinternalGetXvlmMatrix(
              X        = X[rownames(X) != rownames(X)[k], , drop = FALSE],
              formulas = model$formula, 
              parNames = model$model$etaNames
            ),
            coefStart = cf,
            etaStart  = eta[-k, , drop = FALSE] + offset[-k, , drop = FALSE],
            family = model$model,
            priorWeights = pw[-k],
            method = if (is.null(model$call$method)) "IRLS" else model$call$method,
            offset = offset[-k, , drop = FALSE]
          )$beta
        } else {
          kk <- rep(0, length(pw))
          kk[k] <- 1
          
          res[k, ] <- cf - estimatePopsizeFit(
            control = controlMethod(
              silent = TRUE,
              maxiter = maxitNew + 1,
              ...
            ),
            y = y,
            X = model.matrix(model, "vlm"),
            coefStart = cf,
            etaStart  = eta + offset,
            family = model$model,
            priorWeights = pw - kk,
            method = if (is.null(model$call$method)) "IRLS" else model$call$method,
            offset = offset
          )$beta
        }
      }
    }
  }
  
  colnames(res) <- names(cf)
  res
}

#' @method residuals singleRStaticCountData
#' @importFrom stats residuals
#' @rdname regDiagSingleR
#' @exportS3Method
residuals.singleRStaticCountData <- function(object,
                                             type = c(
                                               "pearson",
                                               "pearsonSTD",
                                               "response",
                                               "working",
                                               "deviance",
                                               "all"
                                             ),
                                             ...) {
  type <- match.arg(type)
  res <- object$residuals
  wts <- object$priorWeights
  #mu <- object$fitt.values
  y <- if (is.null(object$y)) stats::model.response(model.frame(object)) else object$y
  
  if (type == "pearsonSTD" && length(object$model$etaNames) > 1)
    stop(paste0("Standardized pearson residuals not yet",
                "implemented for models with multiple linear predictors"))
  
  rs <- switch(
    type,
    working = as.data.frame(
      object$model$funcZ(eta = object$linearPredictors, weight = object$weights, y = y, prior = wts), 
      col.names = paste0("working:", object$model$etaNames)
    ),
    response = res,
    pearson = data.frame(
      "pearson" = res$truncated / sqrt(object$model$variance(eta = object$linearPredictors, type = "trunc"))
    ),
    pearsonSTD = data.frame("pearsonSTD" = res$truncated / sqrt((1 - hatvalues(object)) * object$model$variance(eta = object$linearPredictors, type = "trunc"))),
    deviance = data.frame(
      "deviance" = object$model$devResids(y = y, eta = object$linearPredictors, wt = wts)
    ),
    all = {colnames(res) <- c("truncatedResponse", 
                              "nontruncatedResponse");
      data.frame(
        as.data.frame(object$model$funcZ(eta = object$linearPredictors, weight = object$weights, y = y, prior = wts),
                      col.names = paste0("working:", object$model$etaNames)),
        res,
        "pearson" = as.numeric(res$truncated / sqrt(object$model$variance(eta = object$linearPredictors, type = "trunc"))),
        "pearsonSTD" = if (length(object$model$etaNames) == 1) as.numeric(
          res$truncated / sqrt((1 - hatvalues(object)) * object$model$variance(object$linearPredictors, type = "trunc"))
        ) else NA,
        "deviance" = as.numeric(object$model$devResids(y = y, eta = object$linearPredictors, wt = wts)),
        row.names = rownames(object$linearPredictors)
    )}
  )
  rs
}

#' @importFrom stats cooks.distance
#' @method cooks.distance singleRStaticCountData
#' @rdname regDiagSingleR
#' @exportS3Method 
cooks.distance.singleRStaticCountData <- function(model, ...) {
  if (length(model$model$etaNames) > 1) 
    stop("Cooks distance is only implemented for single parameter families.")
  
  res <- residuals(model, type = "pearsonSTD") ^ 2
  res <- res[, 1]
  
  ht <- hatvalues(model)
  res <- (res * (ht / (length(coef(model)))))
  rownames(res) <- rownames(ht) # fixed
  res
}

# There is no need for documenting the following methods:

#' @method family singleRStaticCountData
#' @importFrom stats family
#' @exportS3Method
family.singleRStaticCountData <- function(object, ...) {
  object$model
}

#' @method print summarysingleRmargin
#' @exportS3Method 
print.summarysingleRmargin <- function(x, ...) {
  cat("Test for Goodness of fit of a regression model:\n",
      "\n", sep = "")
  print(x$Test)
  cat("\n--------------------------------------------------------------",
      "\nCells with fitted frequencies of < 5 have been", x$l5, 
      "\nNames of cells used in calculating test(s) statistic:", names(x$y), 
      "\n", sep = " ")
}

#' @method AIC singleRStaticCountData
#' @importFrom stats AIC
#' @exportS3Method 
AIC.singleRStaticCountData <- function(object, ...) {
  2 * (length(object$coefficients) - object$logL)
}
#' @method BIC singleRStaticCountData
#' @importFrom stats BIC
#' @exportS3Method 
BIC.singleRStaticCountData <- function(object, ...) {
  length(object$coefficients) * log(nobs(object, ...)) - 2 * object$logL
}
#' @method extractAIC singleRStaticCountData
#' @importFrom stats extractAIC
#' @exportS3Method 
extractAIC.singleRStaticCountData <- function(fit, scale, k = 2, ...) {
  -2 * fit$logL + k * length(fit$coefficients)
}
# CODE MODIFIED FROM stats:::logLik.glm
#' @method logLik singleRStaticCountData
#' @importFrom stats logLik
#' @exportS3Method 
logLik.singleRStaticCountData <- function(object, ...) {
  val <- object$logL
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(stats::coef(object))
  class(val) <- "logLik"
  val
}
# CODE MODIFIED FROM stats:::model.frame.glm
#' @method model.frame singleRStaticCountData
#' @importFrom stats glm
#' @importFrom stats model.frame
#' @importFrom stats update
#' @exportS3Method 
model.frame.singleRStaticCountData <- function(formula, ...) {
  dots <- list(...)
  dotargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  
  if (length(dotargs) || is.null(formula$modelFrame)) {
    # fcall <- formula$call
    # fcall$method <- "model.frame"
    # fcall[[1L]] <- quote(stats::glm)
    # fcall[names(nargs)] <- nargs
    # env <- environment(formula$terms)
    # eval(fcall, env)
    # TODO:: low priority add na action and subset here
    combinedFromula <- singleRinternalMergeFormulas(formula$formula)
    if (!is.null(dotargs$data)) {
      jj <- all.vars(combinedFromula)[attr(terms(combinedFromula), "response")]
      if (!(jj %in% colnames(dotargs$data))) {
        combinedFromula <- update(combinedFromula, NULL ~ .)
      }
    }
    stats::model.frame(
      combinedFromula,
      data = if (is.null(dotargs$data)) eval(formula$call$data) else dotargs$data
    )
  }
  else formula$modelFrame
}

#' @method model.matrix singleRStaticCountData
#' @importFrom stats model.matrix
#' @exportS3Method 
model.matrix.singleRStaticCountData <- function(object, type = c("lm", "vlm"), ...) {
  if (missing(type)) type <- "lm"
  
  switch (type,
    lm = {
      X <- model.frame(object);
      X <- model.matrix(object$terms, X)
    },
    vlm = {
      X <- model.frame(object, ...);
      singleRinternalGetXvlmMatrix(
        X = X, 
        formulas = object$formula, 
        parNames = object$model$etaNames
      );
    }
  )
}
#' @method redoPopEstimation singleRStaticCountData
#' @rdname redoPopEstimation
#' @exportS3Method
redoPopEstimation.singleRStaticCountData <- function(object, 
                                                     newdata, 
                                                     cov, 
                                                     weights,
                                                     coef,
                                                     control,
                                                     popVar,
                                                     offset,
                                                     weightsAsCounts,
                                                     ...) {
  ### weightsAsPopCount works
  if (missing(cov)) {
    cov <- vcov
  }
  
  if (is.character(cov)) {
    cov <- get(cov, mode = "function", envir = parent.frame())
  }
  if (is.function(cov)) {
    cov <- cov(object, ...)
  }
  
  if (missing(newdata)) {
    Xvlm <- model.matrix(object, "vlm")
    
    pw <- if (missing(weights))
      object$priorWeights
    else weights
    
    offset <- if (missing(offset))
      object$offset
    else offset
    
    etaNew <- if (missing(coef))
      object$linearPredictors
    else matrix(as.matrix(Xvlm) %*% coef, ncol = length(family(object)$etaNames))
    
    if (missing(control))
      control <- object$populationSize$control
    
    if (is.null(control$bootstrapFitcontrol)) {
      control$bootstrapFitcontrol <- object$control$controlMethod
    }
    
    Y <- object$y
    X <- model.matrix(object)
    MM <- model.frame(object, ...)
    nn <- nobs(object)
    
  } else {
    if (missing(control))
      control <- object$populationSize$control
    
    if (is.null(control$bootstrapFitcontrol)) {
      control$bootstrapFitcontrol <- object$control$controlMethod
    }
    
    MM <- model.frame(object, data = newdata, ...)
    X <- model.matrix(object, data = newdata)
    Xvlm <- model.matrix(object, type = "vlm", data = newdata)
    
    Y <- model.response(MM)
    
    nn <- length(Y)
    
    pw <- if (missing(weights))
      rep(1, nn)
    else weights
    
    offset <- if (missing(offset))
      matrix(0, nrow = length(pw), ncol = length(family(object)$etaNames))
    else offset
    
    coef <- if (missing(coef)) stats::coef(object)
    
    etaNew <- matrix(
      as.matrix(Xvlm) %*% coef, 
      ncol = length(family(object)$etaNames),
      dimnames = list(
        rownames(Xvlm),
        family(object)$etaNames
      )
    ) + offset
  }
  
  singleRcaptureinternalpopulationEstimate(
    y = Y,
    formulas = object$formula,
    X = X,
    grad = object$model$makeMinusLogLike(
      y = Y,
      X = Xvlm,
      weight = pw, 
      deriv = 1
    ),
    hessian = object$model$makeMinusLogLike(
      y = Y,
      X = Xvlm, 
      weight = pw, 
      deriv = 2
    ),
    popVar = if (missing(popVar)) 
      "analytic" 
    else popVar,
    weights = pw,
    eta = etaNew,
    family = family(object),
    beta = if (missing(coef))
      stats::coef(object)
    else coef,
    control = if (missing(control))
      object$populationSize$control
    else control,
    Xvlm = Xvlm,
    W = if (isTRUE(object$call$method == "IRLS") & missing(newdata)) 
      object$weights 
    else 
      family(object)$Wfun(prior = pw, eta = etaNew, y = Y),
    sizeObserved = nn,
    modelFrame = MM,
    cov = cov,
    offset = offset,
    weightsFlag = if (missing(weightsAsCounts)) object$control$controlModel$weightsAsCounts else weightsAsCounts
  )
}
#' @method dfpopsize singleRStaticCountData
#' @rdname regDiagSingleR
#' @exportS3Method 
dfpopsize.singleRStaticCountData <- function(model, dfbeta = NULL, ...) {
  if (isTRUE(model$call$popVar == "bootstrap")) 
    warning("dfpopsize may (in some cases) not work correctly when bootstrap was chosen as population variance estimate.")
  
  dfb <- if (is.null(dfbeta)) dfbeta(model, ...) else dfbeta
  
  X <- model.frame(model, ...)
  
  N <- model$populationSize$pointEstimate
  range <- 1:NROW(dfb)
  
  res <- vector("numeric", length = NROW(dfb))
  pw <- model$priorWeights
  y <- model.response(model.frame(model))
  
  for (k in range) {
    cf <- model$coefficients - dfb[k, ]
    
    if (isTRUE(model$control$controlModel$weightsAsCounts == FALSE)) {
      res[k] <- model$model$pointEst(
        eta = matrix(
          singleRinternalGetXvlmMatrix(
            X = X[rownames(X) != rownames(X)[k], , drop = FALSE], 
            formulas = model$formula, 
            parNames = model$model$etaNames
          ) %*% cf, 
          ncol = length(model$model$etaNames)
        ),
        y = y[-k],
        pw = pw[-k]
      )
    } else {
      # Here additional conditional is not needed since if weights are zero nothing breaks
      kk <- rep(0, length(pw))
      kk[k] <- 1
      res[k] <- model$model$pointEst(
        eta = matrix(
          model.matrix(model, "vlm") %*% cf, 
          ncol = length(model$model$etaNames)
        ),
        y = y,
        pw = pw - kk
      )
    }
  }
  
  N - res
}
#' @method stratifyPopsize singleRStaticCountData
#' @rdname stratifyPopsize
#' @importFrom stats vcov
#' @importFrom stats contrasts
#' @exportS3Method
stratifyPopsize.singleRStaticCountData <- function(object, 
                                                   stratas,
                                                   alpha, 
                                                   cov = NULL,
                                                   ...) {
  ## TODO:: New data doesn't work yet
  
  # if stratas is unspecified get all levels of factors in modelFrame
  if (missing(stratas)) {
    stratas <- names(which(attr(object$terms, "dataClasses") == "factor"))
    stratas <- stratas[stratas %in% attr(object$terms, "term.labels")]
    if (!length(stratas)) {
      stratas <- names(which(attr(object$terms, "dataClasses") == "character"))
      stratas <- stratas[stratas %in% attr(object$terms, "term.labels")]
    }
    if (!length(stratas)) {
      stop("No stratas argument was provided and no factors or character columns are present in model.frame.")
    }
  }
  # If there are no factors or characters and no stratas was provided throw error
  # if significance level is unspecified set it to 5%
  if (missing(alpha)) alpha <- .05
  
  # convert stratas to list for all viable types of specifying the argument
  if (inherits(stratas, "formula")) {
    mf <- model.frame(stratas, model.frame(object))
    mmf <- model.matrix(
      stratas, data = mf, 
      contrasts.arg = lapply(
        subset(mf, select = sapply(mf, is.factor)), # this makes it so that all levels of factors are encoded
        contrasts, contrasts = FALSE
      )
    )
    trm <- attr(mf, "terms")
    stratas <- list()
    for (k in attr(trm, "term.labels")) {
      if (k %in% colnames(mf)) {
        if (is.integer(mf[, k]) | is.character(mf[, k])) {
          for (t in unique(mf[,k])) {
            stratas[[paste0(k, "==", t)]] <- mf[,k] == t
          }
        } else if (is.factor(mf[, k])) {
          for (t in levels(mf[, k])) {
            stratas[[paste0(k, "==", t)]] <- mf[,k] == t
          }
        }
      } else {
        tLevs <- colnames(mmf)[attr(mmf, "assign") == which(attr(trm, "term.labels") == k)]
        for (t in tLevs) {
          stratas[[as.character(t)]] <- mmf[, t] == 1
        }
      }
    }
  } else if (is.list(stratas)) {
    if (!all(sapply(stratas, is.logical)))
      stop("Invalid way of specifying subpopulations in stratas. If stratas argument is a list ")
    
    if (length(stratas[[1]]) != object$sizeObserved) 
      stop("Elements of stratas object should have length equal to number of observed units.")
    
  } else if (is.logical(stratas)) {
    if (length(stratas) != object$sizeObserved) 
      stop("Stratas object should have length equal to number of observed units.")
    
    stratas <- list(strata = stratas)
  } else if (is.character(stratas)) {
    modelFrame <- model.frame(object)
    out <- list()
    for (k in stratas) {
      if (!(k %in% colnames(modelFrame))) 
        stop("Variable specified in stratas is not present in model frame.")
      
      #if (!(is.factor(modelFrame[, k])) & !(is.character(modelFrame[, k]))) 
      #  stop("Variable specified in stratas is not a factor or a character vector.")

      if (is.factor(modelFrame[, k])) {
        # this makes a difference on factor that is not present
        for (t in levels(modelFrame[, k])) {
          out[[paste0(as.character(k), "==", t)]] <- (modelFrame[, k] == t)
        }
      } else {
        for (t in unique(modelFrame[, k])) {
          out[[paste0(as.character(k), "==", t)]] <- (modelFrame[, k] == t)
        }
      }
    }
    stratas <- out
  } else {
    # a formula, a list with logical vectors specifying different sub 
    # populations or a single logical vector or a vector 
    # with names of factor variables.
    errorMessage <- paste0(
      "Invalid way of specifying subpopulations in stratas.\n", 
      "Please provide either:\n",
      "(1) - a list with logical vectors specifying different sub populations\n",
      "(2) - a single logical vector\n",
      "(3) - a formula\n",
      "(4) - a vector with names of variables by which stratas will be created\n"
    )
    stop(errorMessage)
  }
  
  # get necessary model info AFTER possible error in function
  family <- family(object = object)
  priorWeights <- object$priorWeights
  eta <- object$linearPredictors
  Xvlm <- model.matrix(object, "vlm")
  # this is now needed
  y <- if (is.null(object$y)) model.response(model.frame(object)) else object$y
  flagWeighting <- object$control$controlModel$weightsAsCounts
  
  # get covariance matrix
  if (is.function(cov)) cov <- cov(object, ...)
  if (is.null(cov)) cov <- vcov(object, ...)
  
  obs <- vector(mode = "numeric", length = length(stratas))
  est <- vector(mode = "numeric", length = length(stratas))
  stdErr <- vector(mode = "numeric", length = length(stratas))
  cnfStudent <- matrix(nrow = length(stratas), ncol = 2)
  cnfChao <- matrix(nrow = length(stratas), ncol = 2)
  sc <- qnorm(p = 1 - alpha / 2)
  if (length(sc) != length(stratas)) sc <- rep(sc, length.out = length(stratas))
  
  for (k in 1:length(stratas)) {
    cond <- stratas[[k]]
    
    if (isTRUE(object$control$controlModel$weightsAsCounts)) {
      obs[k] <- sum(priorWeights[cond])
    } else {
      obs[k] <- sum(cond)
    }
    
    if (obs[k] > 0) {
      est[k] <- family$pointEst(pw  = priorWeights[cond], 
                                eta = eta[cond, , drop = FALSE],
                                y   = y[cond])
      
      stdErr[k] <- family$popVar(
        pw = priorWeights[cond], 
        eta = eta[cond, , drop = FALSE], 
        cov = cov, 
        Xvlm = subset(Xvlm, subset = rep(cond, length(family$etaNames))),
        y   = y[cond]
      ) ^ .5
      
      cnfStudent[k, ] <- est[k] + c(-sc[k] * stdErr[k], sc[k] * stdErr[k])
      
      G <- exp(sc[k] * sqrt(log(1 + (stdErr[k]^2) / ((est[k] - obs[k]) ^ 2))))
      cnfChao[k, ] <- obs[k] + c((est[k] - obs[k]) / G, (est[k] - obs[k]) * G)
    } else {
      est[k] <- 0
      stdErr[k] <- 0
      cnfStudent[k, ] <- c(0, 0)
      cnfChao[k, ] <- c(0, 0)
    }
  }
  
  result <- data.frame(
    obs, est, 100 * obs / est, stdErr, 
    cnfStudent[, 1], cnfStudent[, 2], 
    cnfChao[, 1], cnfChao[, 2],
    names(stratas), alpha
  )
  
  bounds <- c("LowerBound", "UpperBound")
  
  colnames(result) <- c(
    "Observed", "Estimated", 
    "ObservedPercentage", "StdError", 
    paste0("normal", bounds), 
    paste0("logNormal", bounds), 
    "name", "confLevel"
  )
  
  result
}

#' @method popSizeEst singleRStaticCountData
#' @rdname popSizeEst
#' @exportS3Method
popSizeEst.singleRStaticCountData <- function(object, ...) {
  object$populationSize
}

#' @method print popSizeEstResults
#' @exportS3Method 
print.popSizeEstResults <- function(x, ...) {
  cat("Point estimate: ", x$pointEstimate, 
      "\nVariance: ", x$variance, "\n", 
      (1 - x$control$alpha) * 100, "% confidence intervals:\n", 
      sep = "")
  print(x$confidenceInterval)
  
  invisible(x)
}
#' @method print summarysingleRStaticCountData
#' @importFrom stats printCoefmat
#' @exportS3Method 
print.summarysingleRStaticCountData <- function(x, 
                                                signif.stars = getOption("show.signif.stars"), 
                                                digits = max(3L, getOption("digits") - 3L), 
                                                ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Pearson Residuals:\n")
  print(summary(c(x$residuals[, 1])))
  
  # Start accesing the coefficient information in proper format
  
  cat("\nCoefficients:\n")
  
  cond <- matrix(sapply(
    x$model$etaNames, 
    FUN = function(k) {
      sapply(strsplit(rownames(x$coefficients), 
                      split = ":"), 
             FUN = function(x) x[[length(x)]] == k)
    }
  ), 
  ncol = length(x$model$etaNames), 
  dimnames = list(
    NULL,
    x$model$etaNames
    )
  )
  
  for (k in x$model$etaNames) {
    if (!all(cond[, k] == FALSE)) {
      cat("-----------------------\nFor linear predictors associated with:", k, "\n")
      # Base library has no str_length and no proper sub_string this is a 
      # slightly convoluted way of making do without them
      toPrint <- subset(x$coefficients, cond[,k])
      lengths <- sapply(rownames(toPrint), 
                        function(x) {length(strsplit(x, split = "")[[1]])})
      lK <- length(unlist(strsplit(k, split = "")))
      rownames(toPrint) <- sapply(
        1:nrow(toPrint), 
        function(x) {substr(
          x = rownames(toPrint)[x], 
          start = 1, 
          stop = lengths[[x]] - (1 + lK)
        )}
      )
      #print(toPrint)
      printCoefmat(
        toPrint, digits = digits, 
        signif.stars = signif.stars, 
        signif.legend = if (k == x$model$etaNames[length(x$model$etaNames)]) signif.stars else FALSE, 
        P.values = TRUE, has.Pvalue = TRUE, 
        na.print = "NA", 
        ...
      )
    } else {
      cat("-----------------------\nFor linear predictors associated with:", k, "\n")
      #print(subset(x$coefficients, rowSums(cond) == 0))
      printCoefmat(
        subset(x$coefficients, rowSums(cond) == 0), 
        digits = digits, 
        signif.stars = signif.stars, 
        signif.legend = if (k == x$model$etaNames[length(x$model$etaNames)]) signif.stars else FALSE, 
        P.values = TRUE, has.Pvalue = TRUE, 
        na.print = "NA", 
        ...
      )
    }
  }
  
  cat("\n")
  
  if (!is.null(x$correlation)) {
    corr <- round(x$correlation, digits = 2)
    corr[!lower.tri(corr)] <- ""
    print(corr[-1, -dim(corr)[2]], quote = FALSE)
    cat("\n")
  }
  
  sd <- sqrt(x$populationSize$variance)
  
  if (x$populationSize$control$sd == "normalMVUE") {
    sd <- sd / (sqrt(2 / (x$sizeObserved - 1)) * 
    exp(lgamma(x$sizeObserved / 2) - lgamma((x$sizeObserved - 1) / 2)))
  }
  
  cat(
    "AIC: ", x$aic,
    "\nBIC: ", x$bic,
    "\nResidual deviance: ", x$deviance,
    "\n\nLog-likelihood: ", x$logL, " on ", x$dfResidual, " Degrees of freedom ",
    # optim does not allow for accessing information
    # on number of iterations performed only a number 
    # of calls for gradient and objective function
    if (isTRUE(x$call$method == "optim")) "\nNumber of calls to log-likelihood function: "
    else "\nNumber of iterations: " , x$iter[1], 
    "\n-----------------------",
    "\nPopulation size estimation results: ",
    "\nPoint estimate ", x$populationSize$pointEstimate, 
    "\nObserved proportion: ", round(100 * x$sizeObserved / x$populationSize$pointEstimate, 
                                     digits = 1), 
    "% (N obs = ", x$sizeObserved, ")",
    if (!is.null(x$skew)) "\nBoostrap sample skewness: ", 
    if (!is.null(x$skew)) x$skew, 
    if (!is.null(x$skew)) "\n0 skewness is expected for normally distributed variable\n---",
    if (isTRUE(x$call$popVar == "bootstrap")) "\nBootstrap Std. Error " else "\nStd. Error ", 
    sd, "\n", 
    (1 - x$populationSize$control$alpha) * 100, 
    "% CI for the population size:\n", 
    sep = ""
  )
  
  print(x$populationSize$confidenceInterval)
  
  cat((1 - x$populationSize$control$alpha) * 100, 
      "% CI for the share of observed population:\n", 
      sep = "")
  
  dd <- as.data.frame(x$populationSize$confidenceInterval)
  
  if (ncol(dd) == 1) {
    vctpop <- sort(x$populationSize$confidenceInterval, decreasing = TRUE)
    names(vctpop) <- rev(names(vctpop))
    print(100 * x$sizeObserved / vctpop)
  } else {
    print(data.frame(
      lowerBound = 100 * x$sizeObserved / dd[, 2], 
      upperBound = 100 * x$sizeObserved / dd[, 1],
      row.names = rownames(dd)
    ))
  }
  
  
  invisible(x)
}
#' @method print singleRStaticCountData
#' @exportS3Method 
print.singleRStaticCountData <- function(x, ...) {
  cat("Call: ")
  print(x$call)
  cat("\nCoefficients:\n")
  print(coef(x))
  cat("\nDegrees of Freedom:", x$df.null, 
      "Total (i.e NULL);", x$dfResidual, 
      "Residual")
  
  cat("\nAIC: ", signif(AIC(x)), 
      "\nBIC: ", signif(BIC(x)), 
      "\nResidual deviance: ", signif(x$deviance))
  cat("\n-----------------------------------\n",
      "Population size estimation results:\n",
      sep = "")
  print(x$populationSize)

  invisible(x)
}

#' @method print singleRfamily
#' @exportS3Method 
print.singleRfamily <- function(x, hideFormulas = c(FALSE, TRUE), ...) {
  cat("Family of distributions:", x$family,
      "\nNames of parameters:", x$etaNames, 
      "\nLinks:", attr(x$links, "linkNames"), 
      sep = " ")
  if (!is.null(x$extraInfo)) {
    
    nn <- c("mean", "variance", "popSizeEst", "meanTr", "varianceTr")
    if (any(xx <- is.na(x$extraInfo[nn])))
      warning(cat("The: ", nn[is.na(xx)], sep = " ",
                  "slots are empty for family: ", x$family))
    
    if (!hideFormulas[1])
      cat("\n--\nFormula for mean in this distribution: ", x$extraInfo["mean"],
          "\nFormula for variance in this distribution: ", x$extraInfo["variance"],
          "\nFormula for population size estimation in this distribution: sum(", 
          x$extraInfo["popSizeEst"], ")",
          sep = "")
    
    if (!hideFormulas[2])
      cat("\n\nFormula for mean in truncated distribution: ", x$extraInfo["meanTr"],
          "\nFormula for variance in truncated distribution: ", x$extraInfo["varianceTr"],
          sep = "")
  }
  invisible(x)
}

#' @importFrom stats fitted
#' @method fitted singleRStaticCountData
#' @exportS3Method 
fitted.singleRStaticCountData <- function(object,
                           type = c("truncated", 
                                    "nontruncated",
                                    "all"),
                           ...) {
  if (missing(type)) type <- "truncated"
  switch (type,
    truncated    = object$fittValues$truncated,
    nontruncated = object$fittValues$nontruncated,
    all          = object$fittValues
  )
}

#' @importFrom stats nobs
#' @method nobs singleRStaticCountData
#' @exportS3Method 
nobs.singleRStaticCountData <- function(object, ...) {
  object$sizeObserved
}

#' @importFrom stats df.residual
#' @method df.residual singleRStaticCountData
#' @exportS3Method 
df.residual.singleRStaticCountData <- function(object, ...) {
  object$dfResidual
}

#' @title Generating data in singleRcapture
#' 
#' @description
#' An S3 method for \code{stats::simulate} to handle \code{singleRStaticCountData} and 
#' \code{singleRfamily} classes.
#'
#' @param object an object representing a fitted model.
#' @param nsim a numeric scalar specifying:
#' \itemize{
#'    \item number of response vectors to simulate in \code{simulate.singleRStaticCountData}, defaults to \code{1L}.
#'    \item number of units to draw in \code{simulate.singleRfamily}, defaults to \code{NROW(eta)}.
#' }
#' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).
#' @param truncated logical value indicating whether to sample from truncated or
#' full distribution.
#' @param eta a matrix of linear predictors
#' @param ... additional optional arguments.
#' @return a \code{data.frame} with \code{n} rows and \code{nsim} columns.
#' @seealso [stats::simulate()] [singleRcapture::estimatePopsize()]
#' @examples
#' N <- 10000
#' ###gender <- rbinom(N, 1, 0.2)
#' gender <- rep(0:1, c(8042, 1958))
#' eta <- -1 + 0.5*gender
#' counts <- simulate(ztpoisson(), eta = cbind(eta), seed = 1)
#' df <- data.frame(gender, eta, counts)
#' df2 <- subset(df, counts > 0)
#' ### check coverage with summary
#' mod1 <-  estimatePopsize(
#'   formula       = counts ~ 1 + gender, 
#'   data          = df2, 
#'   model         = ztpoisson, 
#'   controlMethod = list(silent = TRUE)
#' )
#' mod1_sims <- simulate(mod1, nsim=10, seed = 1)
#' colMeans(mod1_sims)
#' mean(df2$counts)
#' @author Maciej Beręsewicz, Piotr Chlebicki
#' @importFrom stats simulate
#' @method simulate singleRStaticCountData
#' @exportS3Method
#' @name simulate
#' @export
simulate.singleRStaticCountData <- function(object, nsim = 1, seed = NULL, ...) {
  n <- nobs(object)
  eta <- object$linearPredictors
  
  # Replicate each row in eta priorWeights number of times
  if (isTRUE(object$control$controlModel$weightsAsCounts)) {
    eta <- matrix(
      unlist(sapply(1:NROW(eta), function(x) {
        rep(eta[x,], object$priorWeights[x])
      })),
      ncol = NCOL(eta),
      byrow = TRUE,
      dimnames = list(
        1:n, colnames(eta)
      )
    )
  }
  
  val <- simulate(
    object    = family(object), 
    seed      = seed, 
    nsim      = n * nsim, 
    eta       = eta, 
    truncated = TRUE
  )

  dim(val) <- c(n, nsim)
  val <- as.data.frame(val)
  names(val) <- paste0("sim_", seq_len(nsim))
  val
}

#' @rdname simulate
#' @importFrom stats simulate
#' @method simulate singleRfamily
#' @exportS3Method
#' @export
simulate.singleRfamily <- function(object, 
                                   nsim,
                                   seed = NULL, 
                                   eta, 
                                   truncated = FALSE, 
                                   ...) {
  if (missing(nsim)) nsim <- NROW(eta)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  object$simulate(nsim, eta, lower = ifelse(truncated, 0, -1))
}
