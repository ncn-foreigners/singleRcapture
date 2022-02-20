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
#' vector of fitted coefficients of regression\cr
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
#' vector of fitted values\cr
#' a list containing population size estimate form PopulationEstimate function\cr
#' model, linear predictors, and qr decomposition of model matrix
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

  prior.weights <- weights
  weights <- 1

  Model_frame <- stats::model.frame(formula, data,  ...)
  Observed <- Model_frame[, 1]
  Variables <- stats::model.matrix(formula, Model_frame, ...)

  if (!family$valideta(start) && !is.null(start)) {
    stop("Invalid start parameter")
  }

  log_like <- family$make_minusloglike(y = Observed,
                                       X = as.matrix(Variables),
                                       weight = prior.weights)
  grad <- family$make_gradient(y = Observed,
                               X = as.matrix(Variables),
                               weight = prior.weights)
  Hess <- family$make_hessian(y= Observed,
                              X = as.matrix(Variables),
                              weight = prior.weights)

  start <- stats::glm(formula = formula,
                      family = stats::poisson(),
                      data = data)$coefficients

  if (family$family %in% c("ztpoisson",
                           "zotpoisson",
                           "chao",
                           "zelterman")) {
    df.reduced <- length(Observed) - dim(Variables)[2]
    if (method == "robust") {
      FITT <- IRLS(Dependent = Observed,
                   Covariates = as.matrix(Variables),
                   eps = 1e-7,
                   family = family,
                   weights = prior.weights,
                   start = start)
      iter <- FITT$iter
      weights <- FITT$Weights
      coefficients <- FITT$Coefficients
    } else if (method == "mle") {
      if (dim(Variables)[2] == 1) {
        FITT <- stats::optimise(f = log_like,
                                interval = c(-5 * start,
                                              5 * start))
        coefficients <- FITT$minimum
        iter <- NULL
      } else {
        FITT <- stats::optim(fn = log_like,
                             par = start,
                             gr = grad)
        coefficients <- FITT$par
        iter <- FITT$counts
      }
    } else {
      stop("Method not implemented")
    }

    coefficients <- as.vector(coefficients)
    names(coefficients) <- colnames(Variables)
    Eta <- as.matrix(Variables) %*% coefficients
    Fitted <- family$linkinv(Eta)
    hess <- Hess(beta = coefficients)
  } else if (family$family %in% c("ztnegbin", "zotnegbin")) {
    df.reduced <- length(Observed) - dim(Variables)[2]
    if (method == "robust") {
      # This is just a starting point there's probably a better choice
      dispersion <- abs(mean(Observed ** 2) - mean(Observed)) / (mean(Observed) ** 2)
      FITT <- IRLS(Dependent = Observed,
                   Covariates = as.matrix(Variables),
                   eps = 1e-7,
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
      weights <- FITT$Weights
      beta <- FITT$Coefficients

      dispersion <- FITT$disp
      coefficients <- c(dispersion, beta)
      Eta <- as.matrix(Variables) %*% beta
      Fitted <- family$linkinv(Eta)
    } else if (method == "mle") {
      dispersion <- abs(mean(Observed ** 2) - mean(Observed)) / (mean(Observed) ** 2)
      start <- c(dispersion, start)
      FITT <- stats::optim(fn = log_like,
                           par = start,
                           gr = grad)
      beta <- FITT$par[-1]
      iter <- FITT$counts
      Eta <- as.matrix(Variables) %*% beta
      Fitted <- family$linkinv(Eta)

      C1 <- FITT$par[1]
      dispersion <- C1
      coefficients <- c(C1, beta)
    } else {
      stop("Method not implemented")
    }

    names(coefficients) <- c("(Dispersion)", colnames(Variables))
    hess <- Hess(ARG = coefficients)
  }

  Eta <- as.numeric(Eta)
  names(Eta) <- rownames(Variables)

  if (sum(diag(-solve(hess)) <= 0) != 0) {
    stop("fitting error analytic hessian is invalid,
         try another model")
  } else {
    Std_err <- sqrt(diag(solve(-hess)))
  }
  Tval <- coefficients / Std_err

  null.deviance <- as.numeric(NULL)
  LOG <- -log_like(coefficients)
  ResRes <- weights * (Observed - Fitted)
  aic <- 2 * length(coefficients) - 2 * LOG
  deviance <- as.numeric(NULL)
  # VGAM::vglm używa przybliżenia rozkładem normalnym standaryzowanym można ewentualnie zamienić
  Pvals <- (stats::pt(q = abs(Tval), df = df.reduced, lower.tail = FALSE) +
            stats::pt(q = -abs(Tval), df = df.reduced, lower.tail = TRUE))
  qr <- qr(Variables)

  Pop <- PopulationEstimate(y = Observed,
                            X = as.data.frame(Variables),
                            Grad = grad,
                            Hess = Hess,
                            Method = pop.ci,
                            weights = prior.weights,
                            Parameter = Fitted,
                            family = family,
                            dispersion = dispersion,
                            beta = coefficients,
                            trcount = trcount)

  Result <- list(y = Observed,
                 formula = formula,
                 call = match.call(),
                 coefficients = coefficients,
                 Standard_Errors = Std_err,
                 Tvalues = Tval,
                 Pvalues = Pvals,
                 null.deviance = null.deviance,
                 model = family,
                 aic = aic,
                 deviance = deviance,
                 prior.weights = prior.weights,
                 weights = weights,
                 residuals = ResRes,
                 LogL = LOG,
                 iter = iter,
                 dispersion = dispersion,
                 df.residual = df.reduced,
                 df.null = length(Observed) - 1,
                 fitted.values = Fitted,
                 Population_Size = Pop,
                 model = Model_frame,
                 linear.predictors = Eta,
                 qr = qr,
                 rank = qr$rank)
  class(Result) <- c("SingleR", "glm", "lm")
  Result
}
