#' @title estimate_popsize
#'
#' @param data Data frame for the regression
#' @param formula Description of model to that is supposed to be fitted
#' @param model Model for regression and population estimate, either a name of the model a function or a family class object returned by appropriate function
#' @param weights Optional object of prior weights used in fitting the model
#' @param subset A logical vector indicating which observations should be used in regression and population size estimation
#' @param na.action TODO
#' @param method Method for fitting values currently supported robust (IRLS) and MaxLikelihood (for now handled by [stats::optim()])
#' @param pop.var A method of constructing confidence interval either analytic or bootstrap
#' where bootstraped confidence interval may be specified in control.pop.var. 
#' There is also the third possible value of noEst which skips the population size estimate alltogether.
#' @param control.method A list indicating parameter to use in population size variance estimation may be constructed with singleRcapture::control.method function
#' @param control.model Not yet implemented
#' @param control.pop.var A list indicating parameter to use in population size variance estimation may be constructed with singleRcapture::control.pop.var function
#' @param modelFrame,x,y logical value indicating whether to return model matrix, dependent vector and model matrix as a part of output
#' @param contrasts Not yet implemented
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
#'  \item{dispersion -- Dispersion parameter if applies.}
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
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats optim
#' @importFrom MASS glm.nb
#' @importFrom stats pnorm
#' @importFrom stats family
#' @export
estimate_popsize <- function(formula,
                             data,
                             model = c("ztpoisson", "ztnegbin",
                                       "ztgeom", "zotpoisson", "ztoipoisson",
                                       "zotnegbin", "ztoinegbin", "zotgeom",
                                       "ztoigeom", "zelterman", "chao"),
                             weights = 1,
                             subset = NULL,
                             na.action = NULL,
                             method = c("mle", "robust"),
                             pop.var = c("analytic",
                                         "bootstrap",
                                         "noEst"),
                             control.method = NULL,
                             control.model = NULL,
                             control.pop.var = NULL,
                             modelFrame = FALSE,
                             x = TRUE,
                             y = TRUE,
                             contrasts = NULL,
                             ...) {
  subset <- parse(text = deparse(substitute(subset)))
  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }
  family <- model
  dispersion <- NULL
  omegaTheta <- NULL
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
                        bootstrapFitcontrol = control.method(epsilon = 1e-3, maxiter = 20, mleMethod = if (grepl(x = family$family, pattern = "negbin")) "Nelder-Mead" else "L-BFGS-B", silent = TRUE, dispGiven = TRUE))
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
  #if (grepl(x = family$family, pattern = "negbin") && control.pop.var$covType == "Fisher") {stop("Obtaining covariance martix from fisher information matrix is not yet aviable for negative binomial models")}
  
  modelFrame <- stats::model.frame(formula, data,  ...)
  variables <- stats::model.matrix(formula, modelFrame, contrasts = contrasts, ...)
  terms <- attr(modelFrame, "terms")
  contrasts <- attr(variables, "contrasts")
  
  subset <- eval(subset, modelFrame)
  if (is.null(subset)) {subset <- TRUE}
  # subset is often in conflict with some packages hence explicit call
  modelFrame <- base::subset(modelFrame, subset = subset)
  
  variables <- base::subset(variables, subset = subset)
  observed <- modelFrame[, 1]
  sizeObserved <- nrow(data) + control.pop.var$trcount

  
  if (!is.null(weights)) {
    weights0 <- prior.weights <- as.numeric(weights)
  } else {
    weights0 <- prior.weights <- 1
  }
  weights <- 1
  
  dataOriginal <- list(y = observed, x = as.matrix(variables), 
                       wp = prior.weights, w = weights0)
  
  if(!all(observed > 0)) {
    stop("Error in function estimate.popsize, data contains zero-counts")
  }

  if (!family$valideta(start) && !is.null(start)) {
    stop("Invalid start parameter")
  }
  
  n1 <- colnames(variables)

  if (family$family %in% c("zotpoisson", "zotnegbin", "zotgeom") &&
      1 %in% unique(observed)) {
    tempdata <- data.frame(observed, prior.weights, variables)
    control.pop.var$trcount <- control.pop.var$trcount + length(observed[observed == 1])
    tempdata <- tempdata[tempdata["observed"] > 1, ]

    observed <- tempdata[, 1]
    prior.weights <- tempdata[, 2]
  
    variables <- data.frame(tempdata[, -c(1, 2)])
    colnames(variables) <- n1
  } else if (family$family == "chao" &&
             !(all(as.numeric(names(table(observed))) %in% c(1, 2)))) {

    tempdata <- data.frame(observed, prior.weights, variables)
    control.pop.var$trcount <- control.pop.var$trcount + length(observed[observed > 2])
    tempdata <- tempdata[tempdata["observed"] == 1 | tempdata["observed"] == 2, ]

    observed <- tempdata[, 1]
    prior.weights <- tempdata[, 2]
    variables <- data.frame(tempdata[, -c(1, 2)])
    colnames(variables) <- n1
  }
  
  if (family$family == "zelterman") {
    # In zelterman model regression is indeed based only on 1 and 2 counts
    # but estimation is based on ALL counts
    tempdata <- data.frame(observed, prior.weights, variables)
    tempdata <- tempdata[tempdata["observed"] == 1 | tempdata["observed"] == 2, ]

    if ((dim(variables)[2] == 1) && ("(Intercept)" %in% colnames(variables))) {
      observed <- tempdata[, 1]
      variables <- data.frame("(Intercept)" = rep(1, length(observed)))
      prior.weights <- tempdata[, 2]
    } else {
      observed <- tempdata[, 1]
      prior.weights <- tempdata[, 2]
      variables <- as.matrix(tempdata[, -c(1, 2)])
    }
    colnames(variables) <- n1
  }

  if(colnames(variables)[1] == "X.Intercept.") {
    colnames(variables)[1] <- "(Intercept)"
  }
  
  dataRegression <- list(y = observed, x = as.matrix(variables, rownames.force = TRUE), 
                         wp = prior.weights, w = weights0)
  if (is.null(control.model$start)) {
    start <- stats::glm.fit(
      x = variables,
      y = observed,
      family = stats::poisson()
    )$coefficients
    if (family$family %in% c("chao", "zelterman")) {
      start[1] <- start[1] + log(1 / 2)
    }
  } else {
    start <- control.method$start
  }
  if (grepl("negbin", family$family)) {
    #start <- MASS::glm.nb(formula = formula, data = data)
    #dispersion <- -start$theta
    #start <- start$coefficients
    if (is.null(control.method$dispersionstart)) {
      if (control.model$alphaIntercept) {
        dispersion <- log(abs(mean(observed ** 2) - mean(observed)) / (mean(observed) ** 2 + .25))
        start <- c(start, "(Intercept):alpha" = dispersion)
      }
    } else {
      #dispersion <- control.method$dispersionstart
      dispersion <- start # TODO: gosh this is terrible pick a better method
      nm <- names(start)
      start <- c(start, dispersion)
      names(start) <- c(nm, paste0(nm, ":alpha"))
    }
  }
  if (grepl("^ztoi", family$family)) {
    #omegaTheta <- (sum(observed == 1) / sum(observed != 1)) / (sum(observed == 1) / sum(observed != 1) - 1)
    if (is.null(control.method$omegaIntercept)) {omegaTheta <- (sum(observed == 1) / sum(observed != 1)) / (sum(observed == 1) / sum(observed != 1) - 1)} else {omegaTheta <- control.method$omegastart / (1 - control.method$omegastart)}
  }
  
  Xvlm <- as.matrix(variables)
  if (family$parNum > 1) {
    if (control.model$alphaIntercept) {
      Xvlm <- matrix(nrow = nrow(as.matrix(variables)) * 2, ncol = ncol(as.matrix(variables)) + 1)
      Xvlm[1:nrow(variables), 1:ncol(variables)] <- as.matrix(variables)
      Xvlm[-(1:nrow(variables)), ] <- Xvlm[, ncol(variables) + 1] <- 0
      Xvlm[-(1:nrow(variables)), ncol(variables) + 1] <- 1
    } else {
      Xvlm <- kronecker(Y = variables, X = diag(1, nrow = 2))
    }
  }
  colnames(Xvlm) <- names(start)
  
  FITT <- estimate_popsize.fit(y = observed,
                               X = variables,
                               Xvlm = Xvlm,
                               family = family,
                               control = control.method,
                               method = method,
                               prior.weights = prior.weights,
                               start = start,
                               dispersion = dispersion,
                               omegaTheta = omegaTheta,
                               ...)
  coefficients <- FITT$beta
  names(coefficients) <- names(start)
  iter <- FITT$iter

  if (grepl("negbin", family$family)) {
    df.reduced <- 2 * length(observed) - length(coefficients)
    logLike <- family$makeMinusLogLike(y = observed, X = Xvlm,
                                       weight = prior.weights)
    grad <- family$makeGradient(y = observed, X = Xvlm, weight = prior.weights, lambdaPredNumber = ncol(variables))
    hessian <- family$makeHessian(y = observed, X = Xvlm,
                                  weight = prior.weights,
                                  lambdaPredNumber = ncol(variables))
    hess <- hessian(coefficients)
    eta <- matrix(Xvlm %*% coefficients, ncol = 2)
    colnames(eta) <- c("log(lambda)", "log(alpha)")
    #colnames(eta) <- paste0(family$link, c("(lambda)", "(alpha)"))
    rownames(eta) <- rownames(variables)
  } else {
    logLike <- family$makeMinusLogLike(y = observed, X = as.matrix(variables),
                                       weight = prior.weights)
    grad <- family$makeGradient(y = observed, X = as.matrix(variables),
                                weight = prior.weights)
    hessian <- family$makeHessian(y = observed, X = as.matrix(variables),
                                  weight = prior.weights)
    hess <- hessian(coefficients)
    
    df.reduced <- length(observed) - length(coefficients)
    
    eta <- as.matrix(variables) %*% coefficients
    colnames(eta) <- paste0(family$link, "(lambda)")
    rownames(eta) <- rownames(variables)
  }
  weights <- FITT$weights
  
  lambda <- family$linkinv(eta[, 1])
  
  if (family$family == "zelterman") {
    lambda <- family$linkinv(as.matrix(dataOriginal$x) %*% coefficients)
  }

  fitt <- data.frame("mu" = family$mu.eta(eta = eta),
                     "link" = family$linkinv(eta[, 1])) # change later link functions in family class to act on matrix eta
  if ((sum(diag(-solve(hess)) <= 0) != 0) && (control.pop.var$covType == "observedInform")) {
    stop("fitting error analytic hessian is invalid, try another model")
  }
  
  null.deviance <- as.numeric(NULL)
  LOG <- -logLike(coefficients)
  resRes <- prior.weights * (observed - fitt)
  if (family$family %in% c("zelterman", "chao")) {resRes <- resRes - 1}
  aic <- 2 * (length(coefficients) - LOG)
  bic <- length(coefficients) * log(length(observed)) - 2 * LOG
  deviance <- sum(family$dev.resids(y = dataRegression$y, 
                                    mu = fitt$link,
                                    disp = dispersion,
                                    wt = prior.weights,
                                    eta = eta) ** 2)

  POP <- signleRcaptureinternalpopulationEstimate(
    y = if ((grepl(x = family$family, pattern = "^zot.*") || family$family == "chao") && (pop.var == "analytic")) dataRegression$y else dataOriginal$y,
    X = if ((grepl(x = family$family, pattern = "^zot.*") || family$family == "chao") && (pop.var == "analytic")) dataRegression$x else dataOriginal$x,
    grad = grad,
    hessian = hessian,
    method = pop.var,
    weights = prior.weights,
    weights0 = weights0,
    lambda = lambda,
    eta = eta,
    family = family,
    dispersion = dispersion,
    omegaTheta = omegaTheta,
    beta = coefficients,
    control = control.pop.var,
    Xvlm = Xvlm
  )
                       
  structure(
    list(
      y = if (isTRUE(y)) {if (family$family %in% c("chao", "zelterman")) dataOriginal$y else dataRegression$y} else NULL,
      X = if (isTRUE(x)) {if (family$family %in% c("chao")) dataRegression$x else dataOriginal$x} else NULL,
      formula = formula,
      call = match.call(),
      coefficients = coefficients,
      control = list(control.model = control.model,
                     control.method = control.method),
      null.deviance = null.deviance,
      model = family,
      aic = aic,
      bic = bic,
      deviance = deviance,
      prior.weights = if (family$family == "zelterman") {weights0} else {prior.weights},
      weights = weights,
      residuals = resRes,
      logL = LOG,
      iter = iter,
      dispersion = dispersion,
      omegaTheta = omegaTheta,
      df.residual = df.reduced,
      df.null = length(observed) - 1,
      fitt.values = fitt,
      populationSize = POP,
      modelFrame = if (isTRUE(modelFrame)) modelFrame else NULL,
      linear.predictors = eta,
      trcount = control.pop.var$trcount,
      sizeObserved = sizeObserved,
      terms = terms,
      contrasts = contrasts,
      na.action = na.action
    ),
    class = c("singleR", "glm", "lm")
  )
}
