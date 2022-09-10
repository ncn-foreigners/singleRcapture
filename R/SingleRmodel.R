#' @title estimate_popsize
#'
#' @param data Data frame or object coercible to data.frame class containing data for the regression and population size estimation.
#' @param formula Formula for the model to be fitted.
#' @param model Model for regression and population estimate. Possible values are \code{"ztpoisson"}, \code{"ztnegbin"}, \code{"ztgeom"}, \code{"zotpoisson"}, \code{"zotnegbin"}, \code{"zotgeom"}, \code{"zelterman"}, or \code{"chao"}.
#' @param weights Optional object of a priori weights used in fitting the model.
#' @param subset A logical vector indicating which observations should be used in regression and population size estimation.
#' @param na.action TODO
#' @param method Method for fitting values currently supported: iteratively reweighted least squares (\code{robust}) and maximum likelihood (\code{mle}).
#' @param pop.var A method of constructing confidence interval either analytic or bootstrap.
#' Bootstrap confidence interval type may be specified in control.pop.var. 
#' There is also the third possible value of noEst which skips the population size estimate alltogether.
#' @param control.method A list indicating parameter to use in population size variance estimation may be constructed with \code{singleRcapture::control.method} function.
#' @param control.model TODO
#' @param control.pop.var A list indicating parameter to use in population size variance estimation may be constructed with \code{singleRcapture::control.pop.var} function.
#' @param modelFrame,x,y Logical value indicating whether to return model matrix, dependent vector and model matrix as a part of output.
#' @param contrasts Not yet implemented.
#' @param ... Arguments to be passed to other methods
#'
#' @return Returns an object of classes inherited from glm containing:\cr
#' @returns
#' \itemize{
#'  \item{y -- Vector of dependent variable if specified at function call.}
#'  \item{X -- Model matrix if specified at function call.}
#'  \item{formula -- Formula provided on call.}
#'  \item{call -- Call matching original input.}
#'  \item{coefficients -- A vector of fitted coefficients of regression.}
#'  \item{control -- A list of control parameters for control.method and control.model, control.pop.var is included in populationSize.}
#'  \item{null.deviance TODO}
#'  \item{model -- Model which estimation of population size and regression was built, object of class family.}
#'  \item{aic -- Akaike information criterion.}
#'  \item{bic -- Bayesian information criterion.}
#'  \item{deviance -- Deviance for the model.}
#'  \item{prior.weights -- Prior weight provided on call.}
#'  \item{weights -- If robust method of estimation was chosen weights returned by IRLS, otherwise same as prior.weights.}
#'  \item{residuals -- Vector of raw residuals.}
#'  \item{logL -- Logarithm likelihood obtained at final iteration.}
#'  \item{iter -- Numbers of iterations performed in fitting or if stats::optim was used number of call to loglikelihhod function.}
#'  \item{df.residuals -- Residual degrees of freedom.}
#'  \item{df.null -- Null degrees of freedom.}
#'  \item{fitt.values -- Data frame of fitted values for both mu (the expected value) and lambda (Poisson parameter).}
#'  \item{populationSize -- a list containing information of population size estimate.}
#'  \item{modelFrame -- Model frame if specified at call.}
#'  \item{linear predictors -- vector of fitted linear predictors.}
#'  \item{trcount -- Number of truncated observations.}
#'  \item{sizeObserved -- Number of observations in original model frame.}
#'  \item{terms -- terms object used.}
#'  \item{contrasts -- contrasts specified in function call.}
#'  \item{na.aciton -- na.action used.}
#' }
#' 
#' @seealso [stats::optim()] [control.method()] [control.pop.var()] [control.model()] [estimate_popsize.fit()]
#'
#' @importFrom stats glm.fit
#' @importFrom stats poisson
#' @importFrom stats binomial
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats optim
#' @importFrom MASS glm.nb
#' @importFrom stats pnorm
#' @importFrom stats family
#' @export
estimate_popsize <- function(formula,
                             data,
                             model = c("ztpoisson", "ztnegbin", "ztgeom", 
                                       "zotpoisson", "ztoipoisson", "ztHurdlepoisson", 
                                       "zotnegbin",  "ztoinegbin", "zotgeom",
                                       "ztoigeom", "zelterman", "chao"),
                             weights = NULL,
                             subset = NULL,
                             na.action = NULL,
                             method = c("mle", "robust"),
                             pop.var = c("analytic",
                                         "bootstrap",
                                         "noEst"),
                             control.method = NULL,
                             control.model = NULL,
                             control.pop.var = NULL,
                             modelFrame = TRUE,
                             x = TRUE,
                             y = TRUE,
                             contrasts = NULL,
                             ...) {
  if (missing(method)) method <- "mle"
  if (missing(pop.var)) pop.var <- "analytic"
  subset <- parse(text = deparse(substitute(subset)))
  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }
  family <- model
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  # adding control parameters that may possibly be missing
  # since passing simple lists as control arguments is allowed
  m1 <- control.pop.var
  m1 <- m1[sapply(m1, is.null) == FALSE]
  m2 <- control.pop.var(fittingMethod = match.arg(method), 
                        bootstrapFitcontrol = control.method(epsilon = 1e-3, maxiter = 20, mleMethod = if (grepl(x = family$family, pattern = "negbin") || grepl(x = family$family, pattern = "^ztoi") || grepl(x = family$family, pattern = "^oizt")) "Nelder-Mead" else "L-BFGS-B", silent = TRUE))
  m2 <- m2[names(m2) %in% names(m1) == FALSE]
  control.pop.var <- append(m1, m2)
  m1 <- control.method
  m2 <- control.method(mleMethod = if (grepl(x = family$family, pattern = "negbin") || grepl(x = family$family, pattern = "^ztoi")) "Nelder-Mead" else "L-BFGS-B")
  m2 <- m2[names(m2) %in% names(m1) == FALSE]
  control.method <- append(m1, m2)
  m1 <- control.model
  m2 <- control.model()
  m2 <- m2[names(m2) %in% names(m1) == FALSE]
  control.model <- append(m1, m2)
  
  modelFrame1 <- stats::model.frame(formula, data,  ...)
  variables <- stats::model.matrix(formula, modelFrame1, contrasts = contrasts, ...)
  terms <- attr(modelFrame1, "terms")
  contrasts <- attr(variables, "contrasts")
  
  subset <- eval(subset, modelFrame1)
  if (is.null(subset)) {subset <- TRUE}
  # subset is often in conflict with some common packages hence explicit call
  modelFrame1 <- base::subset(modelFrame1, subset = subset)
  variables <- base::subset(variables, subset = subset)
  observed <- modelFrame1[, 1]
  sizeObserved <- nrow(data) + control.pop.var$trcount

  
  if (!is.null(weights)) {
    prior.weights <- as.numeric(weights)
  } else {
    prior.weights <- rep(1, nrow(modelFrame1))
  }
  weights <- 1
  
  if(!all(observed > 0)) {
    stop("Error in function estimate.popsize, data contains zero-counts")
  }

  # if (!family$valideta(start) && !is.null(start)) {
  #   stop("Invalid start parameter")
  # }
  
  formulas <- list(formula)
  if ("alpha" %in% family$etaNames) {
    formulas <- append(x = formulas, control.model$alphaFormula)
  }
  if ("omega" %in% family$etaNames) {
    formulas <- append(x = formulas, control.model$omegaFormula)
  }
  if ("pi" %in% family$etaNames) {
    formulas <- append(x = formulas, control.model$piFormula)
  }

  wch <- singleRcaptureinternalDataCleanupSpecialCases(family = family, observed = observed, pop.var = pop.var)

  control.pop.var$trcount <- control.pop.var$trcount + wch$trr
  
  # TODO::
  ## move this to family functions
  if (is.null(control.method$start)) {
    start <- stats::glm.fit(
      x = variables[wch$reg, ],
      y = observed[wch$reg],
      family = stats::poisson(),
      weights = prior.weights[wch$reg]
    )$coefficients
    if (isTRUE(control.method$useZtpoissonAsStart)) {
      start <- estimate_popsize.fit(
        y = observed[wch$reg],
        X = variables[wch$reg, ],
        family = ztpoisson(),
        start = start,
        hwm = ncol(variables),
        control = control.method(),
        method = method,
        prior.weights = prior.weights
      )$beta
    }
    if (family$family %in% c("chao", "zelterman")) {
      start[1] <- start[1] + log(1 / 2)
    }
  } else {
    start <- control.method$start
  }
  
  if ("omega" %in% family$etaNames) {
    if (is.null(control.method$omegaStart)) {
      if (control.model$omegaFormula == ~ 1) {
        omg <- (length(observed[wch$reg]) - sum(observed == 1)) / (sum(table(observed[wch$reg]) * as.numeric(names(table(observed[wch$reg])))) - length(observed[wch$reg]))
        #start <- c(start, log(omg / (1 - omg)))
        start <- c(start, log(omg))
      } else {
        start <- c(start, stats::glm.fit(
          x = model.matrix(control.model$omegaFormula, modelFrame1[wch$reg, attr(modelFrame1, "names")[-1]]),
          y = as.numeric(observed[wch$reg] == 1),
          family = stats::binomial()
        )$coefficients)
      }
    } else {
      start <- c(start, control.method$omegaStart)
    }
  }
  Xvlm <- singleRinternalGetXvlmMatrix(X = modelFrame1[wch$reg, attr(modelFrame1, "names")[-1]],
                                       nPar = family$parNum, formulas = formulas, parNames = family$etaNames)
  hwm <- Xvlm[[2]]
  Xvlm <- Xvlm[[1]]
  if ("alpha" %in% family$etaNames) {
    if (is.null(control.method$alphaStart)) {
      if (control.model$alphaFormula == ~ 1) {
        start <- c(start, log(abs(mean(observed[wch$reg] ** 2) - mean(observed[wch$reg])) / (mean(observed[wch$reg]) ** 2 + .25)))
      } else {
        cc <- colnames(Xvlm)
        cc <- cc[grepl(x = cc, pattern = "alpha$")]
        cc <- unlist(strsplit(x = cc, ":"))
        cc <- cc[cc != "alpha"]
        start <- c(start, start[cc])  # TODO: gosh this is terrible pick a better method
      }
    } else {
      start <- c(start, control.method$alphaStart)
    }
  }
  if ("pi" %in% family$etaNames) {
    if (is.null(control.method$piStart)) {
      # maybe there is a less complicated way
      cc <- colnames(Xvlm)
      cc <- cc[grepl(x = cc, pattern = "pi$")]
      cc <- unlist(strsplit(x = cc, ":"))
      cc <- cc[cc != "pi"]
      start <- c(start, start[cc])
    } else {
      start <- c(start, control.method$piStart)
    }
  }
  names(start) <- colnames(Xvlm)
  FITT <- estimate_popsize.fit(
    y = observed[wch$reg],
    X = Xvlm,
    family = family,
    control = control.method,
    method = method,
    prior.weights = prior.weights[wch$reg],
    start = start,
    hwm = hwm,
    ...
  )
  coefficients <- FITT$beta
  names(coefficients) <- names(start)
  iter <- FITT$iter
  df.reduced <- nrow(Xvlm) - length(coefficients)
  logLike <- family$makeMinusLogLike(y = observed[wch$reg], X = Xvlm,
                                     weight = prior.weights[wch$reg])
  grad <- family$makeGradient(y = observed[wch$reg], X = Xvlm, weight = prior.weights[wch$reg])
  hessian <- family$makeHessian(y = observed[wch$reg], X = Xvlm,
                                weight = prior.weights[wch$reg],
                                lambdaPredNumber = hwm[1])

  hess <- hessian(coefficients)
  eta <- matrix(as.matrix(Xvlm) %*% coefficients, ncol = family$parNum)
  colnames(eta) <- family$etaNames
  rownames(eta) <- rownames(variables[wch$reg])
  weights <- FITT$weights
  
  if (family$family == "zelterman") {
    eta <- matrix(as.matrix(variables) %*% coefficients, ncol = 1)
    colnames(eta) <- family$etaNames
    rownames(eta) <- rownames(variables)
  }

  fitt <- data.frame(family$mu.eta(eta = eta),
                     family$mu.eta(eta = eta, type = "nontrunc")) # change later link functions in family class to act on matrix eta
  colnames(fitt) <- c("mu", "link")
  if (control.pop.var$covType == "observedInform") {
    if ((sum(diag(-solve(hess)) <= 0) != 0)) {
      stop("Fitting error observed information matrix obtained from analytic hessian is invalid i.e not positive defined, try another model.")
    }
  }
  
  null.deviance <- as.numeric(NULL)
  LOG <- -logLike(coefficients)
  resRes <- prior.weights * (observed[wch$reg] - fitt)
  if (family$family %in% c("zelterman", "chao")) {resRes <- resRes - 1}
  aic <- 2 * (length(coefficients) - LOG)
  bic <- length(coefficients) * log(length(observed[wch$reg])) - 2 * LOG
  deviance <- sum(family$dev.resids(y = observed[wch$reg], 
                                    wt = prior.weights[wch$reg],
                                    eta = if (family$family == "zelterman") eta[wch$reg] else eta) ** 2)
  POP <- singleRcaptureinternalpopulationEstimate(
    #y = if ((grepl(x = family$family, pattern = "^zot.*") || family$family == "chao") && (pop.var == "analytic")) observed[wch$reg] else observed,
    y = observed[wch$est],
    #X = if ((grepl(x = family$family, pattern = "^zot.*") || family$family == "chao") && (pop.var == "analytic")) Xvlm else variables,
    X = variables[wch$est,],
    grad = grad,
    hessian = hessian,
    pop.var = pop.var,
    weights = prior.weights[wch$est],
    eta = eta,
    family = family,
    beta = coefficients,
    control = control.pop.var,
    hwm = hwm,
    Xvlm = Xvlm,
    W = if (method == "robust") weights else family$Wfun(prior = prior.weights, eta = eta)
  )
  structure(
    list(
      y = observed,
      X = variables,
      formula = formulas,
      call = match.call(),
      coefficients = coefficients,
      control = list(control.model = control.model,
                     control.method = control.method),
      null.deviance = null.deviance,
      model = family,
      aic = aic,
      bic = bic,
      deviance = deviance,
      prior.weights = prior.weights,
      weights = weights,
      residuals = resRes,
      logL = LOG,
      iter = iter,
      df.residual = df.reduced,
      df.null = length(observed) - 1,
      fitt.values = fitt,
      populationSize = POP,
      modelFrame = if (isTRUE(modelFrame)) modelFrame1 else NULL,
      linear.predictors = eta,
      trcount = control.pop.var$trcount,
      sizeObserved = sizeObserved,
      terms = terms,
      contrasts = contrasts,
      na.action = na.action,
      which = wch
    ),
    class = c("singleR", "glm", "lm")
  )
}
