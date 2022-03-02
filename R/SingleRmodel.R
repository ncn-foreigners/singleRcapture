#' @title estimate_popsize
#'
#' @param data Data frame for the regression
#' @param formula Description of model to that is supposed to be fitted
#' @param model Model for regression and population estimate
#' @param weights Optional object of a-priori weights used in fitting the model
#' @param subset As in lm/glm
#' @param na.action ||
#' @param method Method for fitting values currently supported IRLS and MaxLikelihood
#' @param pop.ci A method of constructing confidence interval either analytic or bootstrap
#' where bootstraped confidence interval may either be based on 2.5%-97.5%
#' percientiles ("bootstrapPerc") or studentized CI ("bootstrapSD")
#' @param trcount Optional parameter for Zero-one truncated models, if population estimate
#' is for one inflated model then it specifies one counts and includes them in
#' final population estimate both point and interval, and for zeltermann/chao
#' estimator where it specifies counts of not used in estimate
#' @param control.method ||
#' @param control.model ||
#' @param control.pop.var ||
#' @param model.matrix If true returns model matrix as a part of returned object
#' @param x ||
#' @param y ||
#' @param contrasts ||
#' @param ... Arguments to be passed to other methods like subset in stats::model.frame
#'
#' @return Returns an object of classes inherited from glm containing:\cr
#' @returns
#' y argument used in regression\cr
#' formula provided on call\cr
#' Call on which model is based\cr
#' vector of fitt coefficients of regression\cr
#' vector of Standard errors of coeffcients\cr
#' vectors of test and p-values \cr
#' model used in regression\cr
#' akane information criterion\cr
#' vector of prior weights assingned on call\cr
#' vector of weights if IRLS was the fit method then its a vector what was returned by IRLS\cr
#' vector of working residuals\cr
#' Log-likelihood for the model\cr
#' number of itteration performed in fitting\cr
#' total and reduced degrees of freedom\cr
#' vector of fitt values\cr
#' a list containing population size estimate form PopulationEstimate function\cr
#' model, linear predictors, and qr decomposition of model matrix
#' @importFrom stats glm.fit
#' @importFrom stats poisson
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats optim
#' @importFrom MASS glm.nb
#' @importFrom stats pt
#' @export
estimate_popsize <- function(formula,
                             data,
                             model = c("ztpoisson", "ztnegbin",
                                       "zotpoisson", "zotnegbin",
                                       "zelterman", "chao"),
                             weights = 1,
                             subset = NULL,
                             na.action = NULL,
                             trcount = 0,
                             method = c("mle", "robust"),
                             pop.ci = c("analytic",
                                        "bootstrapPerc",
                                        "bootstrapSD"),
                             control.method = NULL,
                             control.model = NULL,
                             control.pop.var = NULL,
                             model.matrix = TRUE,
                             x = FALSE,
                             y = FALSE,
                             contrasts = NULL,
                             ...) {
  family <- model
  dispersion <- NULL
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  model_frame <- stats::model.frame(formula, data,  ...)
  observed <- model_frame[, 1]
  variables <- stats::model.matrix(formula, model_frame, ...)

  if (!family$valideta(start) && !is.null(start)) {
    stop("Invalid start parameter")
  }

  if (!is.null(weights)) {
    prior.weights <- as.numeric(weights)
  } else {
    prior.weights <- 1 / length(observed)
  }
  weights <- 1

  log_like <- family$make_minusloglike(y = observed,
                                       X = as.matrix(variables),
                                       weight = prior.weights)
  grad <- family$make_gradient(y = observed,
                               X = as.matrix(variables),
                               weight = prior.weights)
  hessian <- family$make_hessian(y= observed,
                              X = as.matrix(variables),
                              weight = prior.weights)

  if (family$family %in% c("ztpoisson", "zotpoisson",
                           "chao", "zelterman")) {
    start <- stats::glm.fit(x = variables,
                            y = observed,
                            family = stats::poisson())$coefficients
    if (family$family %in% c("chao", "zelterman")) {
      start[1] <- start[1] + log(1 / 2)
    }
  } else {
    start <- MASS::glm.nb(formula = formula, data = data)
    dispersion <- -start$theta
    start <- start$coefficients
  }

  if (family$family %in% c("ztpoisson",
                           "zotpoisson",
                           "chao")) {
    df.reduced <- length(observed) - dim(variables)[2]
    if (method == "robust") {
      FITT <- IRLS(dependent = observed,
                   covariates = as.matrix(variables),
                   eps = .Machine$double.eps,
                   family = family,
                   weights = prior.weights,
                   start = start)
      iter <- FITT$iter
      weights <- FITT$weights
      coefficients <- FITT$coefficients
    } else if (method == "mle") {
      if (dim(variables)[2] == 1) {
        FITT <- stats::optim(par = start,
                             lower = start - abs(start) * 10,
                             upper = start + abs(start) * 10,
                             fn = log_like,
                             gr = function(x) -grad(x),
                             method = "Brent",
                             control = list(reltol = .Machine$double.eps))
        coefficients <- FITT$par
        iter <- FITT$counts
      } else {
        FITT <- stats::optim(fn = log_like,
                             par = start,
                             gr = function(x) -grad(x),
                             method = "L-BFGS-B",
                             control = list(factr = .Machine$double.eps))
        coefficients <- FITT$par
        iter <- FITT$counts
      }
    } else {
      stop("Method not implemented")
    }

    coefficients <- as.vector(coefficients)
    names(coefficients) <- colnames(variables)
    eta <- as.matrix(variables) %*% coefficients
    fitt <- family$linkinv(eta)
    hess <- hessian(beta = coefficients)
  } else if (family$family %in% c("ztnegbin", "zotnegbin")) {
    df.reduced <- length(observed) - dim(variables)[2]
    if (method == "robust") {
      FITT <- IRLS(dependent = observed,
                   covariates = as.matrix(variables),
                   eps = .Machine$double.eps,
                   disp = dispersion,
                   family = family,
                   weights = prior.weights,
                   start = start)
      if (is.null(FITT)) {
        return(NULL)
        stop("fitting error try another model
             (negative binomial models are highly volitile)")
      }
      iter <- FITT$iter
      weights <- FITT$W
      beta <- FITT$coefficients

      dispersion <- FITT$disp
      coefficients <- c(dispersion, beta)
      eta <- as.matrix(variables) %*% beta
      fitt <- family$linkinv(eta)
    } else if (method == "mle") {
      start <- c(dispersion, start)
      FITT <- stats::optim(fn = log_like,
                           par = start,
                           gr = function(x) -grad(x),
                           method = "L-BFGS-B",
                           control = list(factr = .Machine$double.eps))
      beta <- FITT$par[-1]
      iter <- FITT$counts
      eta <- as.matrix(variables) %*% beta
      fitt <- family$linkinv(eta)

      dispersion <- FITT$par[1]
      coefficients <- c(dispersion, beta)
    } else {
      stop("Method not implemented")
    }
    names(coefficients) <- c("(Dispersion)", colnames(variables))
    hess <- hessian(arg = coefficients)
  } else if (family$family == "zelterman") {
    # In zelterman model regression is indeed based only on 1 and 2 counts
    # but estimation is based on ALL counts
    name1 <- colnames(model_frame)[1]
    tempdata <- model_frame[model_frame[name1] == 1 | model_frame[name1] == 2, ]
    if (is.vector(tempdata)) {
      tempobserved <- tempdata
      tempvariables <- data.frame("(Intercept)" = rep(1, length(tempdata)))
    } else {
      tempobserved <- tempdata[, 1]
      tempvariables <- stats::model.matrix(formula, tempdata, ...)
    }

    log_like <- family$make_minusloglike(y = tempobserved,
                                         X = as.matrix(tempvariables),
                                         weight = prior.weights)

    grad <- family$make_gradient(y = tempobserved,
                                 X = as.matrix(tempvariables),
                                 weight = prior.weights)
    hessian <- family$make_hessian(y= tempobserved,
                                   X = as.matrix(tempvariables),
                                   weight = prior.weights)

    df.reduced <- length(tempobserved) - dim(tempvariables)[2]
    if (method == "robust") {
      FITT <- IRLS(dependent = tempobserved,
                   covariates = as.matrix(tempvariables),
                   eps = .Machine$double.eps,
                   family = family,
                   weights = prior.weights,
                   start = start)
      iter <- FITT$iter
      weights <- FITT$weights
      coefficients <- FITT$coefficients
    } else if (method == "mle") {
      if (dim(variables)[2] == 1) {
        FITT <- stats::optim(par = start,
                             lower = start - abs(start) * 10,
                             upper = start + abs(start) * 10,
                             fn = log_like,
                             gr = function(x) -grad(x),
                             method = "Brent",
                             control = list(reltol = .Machine$double.eps))
        coefficients <- FITT$par
        iter <- FITT$counts
      } else {
        FITT <- stats::optim(fn = log_like,
                             par = as.vector(start),
                             gr = function(x) -grad(x),
                             method = "L-BFGS-B",
                             control = list(factr = .Machine$double.eps))
        coefficients <- FITT$par
        iter <- FITT$counts
      }
    } else {
      stop("Method not implemented")
    }

    coefficients <- as.vector(coefficients)
    names(coefficients) <- colnames(variables)
    eta <- as.matrix(variables) %*% coefficients
    fitt <- family$linkinv(eta)
    hess <- hessian(beta = coefficients)
  }

  eta <- as.numeric(eta)
  names(eta) <- rownames(variables)

  if (sum(diag(-solve(hess)) <= 0) != 0) {
    stop("fitting error analytic hessian is invalid,
         try another model")
  } else {
    stdErr <- sqrt(diag(solve(-hess)))
  }
  wVal <- coefficients / stdErr

  null.deviance <- as.numeric(NULL)
  LOG <- -log_like(coefficients)
  resRes <- weights * (observed - fitt)
  aic <- 2 * length(coefficients) - 2 * LOG
  deviance <- as.numeric(NULL)
  # VGAM::vglm używa przybliżenia rozkładem normalnym standaryzowanym można ewentualnie zamienić
  pVals <- (stats::pt(q =  abs(wVal), df = df.reduced, lower.tail = FALSE) +
            stats::pt(q = -abs(wVal), df = df.reduced, lower.tail = TRUE))
  qr <- qr(variables)

  POP <- populationEstimate(y = observed,
                            X = as.data.frame(variables),
                            grad = grad,
                            hessian = hessian,
                            method = pop.ci,
                            weights = prior.weights,
                            parameter = fitt,
                            family = family,
                            dispersion = dispersion,
                            beta = coefficients,
                            trcount = trcount)

  result <- list(y = observed,
                 formula = formula,
                 call = match.call(),
                 coefficients = coefficients,
                 standard_errors = stdErr,
                 wValues = wVal,
                 pValues = pVals,
                 null.deviance = null.deviance,
                 model = family,
                 aic = aic,
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
                 qr = qr,
                 rank = qr$rank)
  class(result) <- c("singleR", "glm", "lm")
  result
}
