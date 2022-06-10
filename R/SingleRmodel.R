#' @title estimate_popsize
#'
#' @param data Data frame for the regression
#' @param formula Description of model to that is supposed to be fitted
#' @param model Model for regression and population estimate
#' @param weights Optional object of a-priori weights used in fitting the model
#' @param subset Same as in glm
#' @param na.action TODO
#' @param method Method for fitting values currently supported IRLS and MaxLikelihood
#' @param pop.var A method of constructing confidence interval either analytic or bootstrap
#' where bootstraped confidence interval may either be based on 2.5%-97.5%
#' percientiles ("bootstrapPerc") or studentized CI ("bootstrapSD")
#' @param control.method ||
#' @param control.model ||
#' @param control.pop.var A list indicating parameter to use in population size variance estimation
#' @param model.matrix If true returns model matrix as a part of returned object
#' @param x manually provided model matrix
#' @param y manually provided vector for dependent variable
#' @param contrasts ||
#' @param ... Arguments to be passed to other methods
#'
#' @return Returns an object of classes inherited from glm containing:\cr
#' @returns
#' \itemize{
#'  \item{y argument of observations}
#'  \item{formula provided on call}
#'  \item{call on which model is based}
#'  \item{coefficients a vector of fitted coefficients of regression}
#'  \item{standard_errors SE of fitted coefficients of regression}
#'  \item{wValues test values for Wald test}
#'  \item{pValues P values for Wald test}
#'  \item{null.deviance TODO}
#'  \item{model on which estimation and regression was built}
#'  \item{aic akaike information criterion}
#'  \item{bic bayesian information criterion}
#'  \item{deviance TODO}
#'  \item{prior.weight weight provided on call}
#'  \item{weights if robust method of estimation was chosen, weights returned by IRLS}
#'  \item{residuals working residuals}
#'  \item{logL logarithm likelihood obtained}
#'  \item{iter numbers of iterations performed in fitting}
#'  \item{dispersion dispersion parameter if applies}
#'  \item{df.residuals residual degrees of freedom}
#'  \item{df.null null degrees of freedom}
#'  \item{fitt.values vector of fitted values}
#'  \item{populationSize a list containing information of population size estimate}
#'  \item{linear predictors vector of linear predictors}
#'  \item{qr qr decomposision of model matrix, might be removed later since it's
#'  not used anywhere in package}
#' }
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
                                       "zotpoisson", "zotnegbin",
                                       "zelterman", "chao",
                                       "ztgeom", "zotgeom"),
                             weights = 1,
                             subset = NULL,
                             na.action = NULL,
                             method = c("mle", "robust"),
                             pop.var = c("analytic",
                                         "bootstrap"),
                             control.method = NULL,
                             control.model = NULL,
                             control.pop.var = NULL,
                             model.matrix = TRUE,
                             x = FALSE,
                             y = FALSE,
                             contrasts = NULL,
                             ...) {
  subset <- parse(text = deparse(substitute(subset)))
  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }
  # adding control parameters that may possibly be missing
  m1 <- control.pop.var
  m2 <- control.pop.var()
  m2 <- m2[names(m2) %in% names(m1) == FALSE]
  control.pop.var <- append(m1, m2)
  m1 <- control.method
  m2 <- control.method()
  m2 <- m2[names(m2) %in% names(m1) == FALSE]
  control.method <- append(m1, m2)
  m1 <- control.model
  m2 <- control.model()
  m2 <- m2[names(m2) %in% names(m1) == FALSE]
  control.model <- append(m1, m2)
  typefitt <- control.model$typefitted
  typefitt <- ifelse(is.null(typefitt), "link", typefitt)
  family <- model
  dispersion <- NULL
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  model_frame <- stats::model.frame(formula, data,  ...)
  variables <- stats::model.matrix(formula, model_frame, ...)
  subset <- eval(subset, model_frame)
  if (is.null(subset)) {subset <- TRUE}
  # subset is often in conflict with some packages hence explicit call
  model_frame <- base::subset(model_frame, subset = subset)
  if (isFALSE(x)) {
    variables <- base::subset(variables, subset = subset)
  } else {
    variables <- x
  }
  if (isFALSE(y)) {
    observed <- model_frame[, 1]
  } else {
    observed <- y
  }
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
  
  dataRegression <- list(y = observed, x = as.matrix(variables), 
                         wp = prior.weights, w = weights0)

  if (family$family %in% c("ztpoisson", "zotpoisson",
                           "chao", "zelterman",
                           "ztgeom", "zotgeom")) {
    start <- stats::glm.fit(x = variables,
                            y = observed,
                            family = stats::poisson())$coefficients
    if (family$family %in% c("chao", "zelterman")) {
      start[1] <- start[1] + log(1 / 2)
    }
  } else {
    #start <- MASS::glm.nb(formula = formula, data = data)
    #dispersion <- -start$theta
    #start <- start$coefficients
    start <- stats::glm.fit(x = variables,
                            y = observed,
                            family = stats::poisson())$coefficients
    dispersion <- log(abs(mean(observed ** 2) - mean(observed)) / (mean(observed) ** 2))
  }

  FITT <- estimate_popsize.fit(y = observed,
                               X = variables,
                               family = family,
                               control = control.method,
                               method,
                               prior.weights = prior.weights,
                               start = c(dispersion, start),
                               dispersion = dispersion,
                               ...)
  coefficients <- FITT$beta

  log_like <- FITT$ll
  grad <- FITT$grad
  hessian <- FITT$hessian

  hess <- FITT$hess
  iter <- FITT$iter
  df.reduced <- FITT$degf
  weights <- FITT$weights

  if (is.null(dispersion)) {
    eta <- as.matrix(variables) %*% coefficients
  } else {
    eta <- as.matrix(variables) %*% coefficients[-1]
  }
  
  lambda <- family$linkinv(eta)
  
  if (family$family == "zelterman") {
    lambda <- family$linkinv(as.matrix(dataOriginal$x) %*% coefficients)
  }

  # Here do fitt a list or data.frame
  # fitt <- data.fram("mu" = family$mu.eta(eta, disp = dispersion))
  fitt <- data.frame("mu" = family$mu.eta(eta, disp = dispersion),
                     "link" = family$linkinv(eta))
  

  if (!is.null(dispersion)) {
    dispersion <- coefficients[1]
  }

  if (sum(diag(-solve(hess)) <= 0) != 0) {
    stop("fitting error analytic hessian is invalid,
         try another model")
  } else {
    stdErr <- sqrt(diag(solve(-hess)))
  }
  wVal <- coefficients / stdErr

  null.deviance <- as.numeric(NULL)
  LOG <- -log_like(coefficients)
  resRes <- prior.weights * (observed - fitt)
  if (family$family %in% c("zelterman", "chao")) {resRes <- resRes - 1}
  aic <- 2 * (length(coefficients) - LOG)
  bic <- length(coefficients) * log(length(observed)) - 2 * LOG
  deviance <- sum(family$dev.resids(y = dataRegression$y, 
                                    mu = fitt$link,
                                    disp = dispersion,
                                    wt = prior.weights) ** 2)
  # In wald W-values have N(0,1) distributions (asymptotically) pnorm is symmetric wrt 0
  pVals <- 2 * stats::pnorm(q =  abs(wVal), lower.tail = FALSE)

  POP <- populationEstimate(y = if ((grepl(x = family$family, pattern = "^zot.*") || family$family == "chao") && (pop.var == "analytic")) dataRegression$y else dataOriginal$y,
                            X = if ((grepl(x = family$family, pattern = "^zot.*") || family$family == "chao") && (pop.var == "analytic")) dataRegression$x else dataOriginal$x,
                            grad = grad,
                            hessian = hessian,
                            method = pop.var,
                            weights = prior.weights,
                            weights0 = weights0,
                            lambda = lambda,
                            family = family,
                            dispersion = dispersion,
                            beta = coefficients,
                            control = control.pop.var)
  structure(
    list(
      y = if (family$family %in% c("chao", "zelterman")) dataOriginal$y else dataRegression$y,
      X = if (isTRUE(model.matrix)) {if (family$family %in% c("chao", "zelterman")) dataRegression$x else dataOriginal$x} else NULL,
      formula = formula,
      call = match.call(),
      coefficients = coefficients,
      standard_errors = stdErr,
      control = list(control.model = control.model,
                     control.method = control.method,
                     control.pop.var = control.pop.var),
      wValues = wVal,
      pValues = pVals,
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
      dispersion = dispersion,
      df.residual = df.reduced,
      df.null = length(observed) - 1,
      fitt.values = fitt,
      populationSize = POP,
      model = model_frame,
      linear.predictors = eta,
      trcount = control.pop.var$trcount
    ),
    class = c("singleR", "glm", "lm")
  )
}
