#' @title SingleRmodel
#'
#' @param data Data frame for the regression
#' @param formula Description of model to that is supposed to be fitted
#' @param start Optional argument specifying initial values of regression parameters
#' @param Fit.Method Method for fitting values currently supported IRLS and MaxLikelihood
#' @param family Family of distributions used in regression
#' @param CI A method of constructing confidence interval either analytic or bootstrap
#' where bootstraped confidence interval may either be based on 2.5%-97.5%
#' percientiles ("bootstrapPerc") or studentized CI ("bootstrapSD")
#' @param dispersion optional parameter if family exhibits over/underdispersion
#' and it has already known/been estimated providing it should improve model fit considerably
#' @param prior.weights Optional object of a-priori weights used in fitting the model
#' @param ... Arguments to be passed to other methods like subset in stats::model.frame
#'
#' @return Returns an object of classes inherited from glm containing:\cr
#' @returns
#' y argument used in regression\cr
#' formula provided on call\cr
#' Call on which model is based\cr
#' vector of fitted coefficients of regression\cr
#' vector of Standard errors of coeffcients\cr
#' vectors of test and p - values \cr
#' family used in regression\cr
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
Single_Rmodel <- function(data,
                          formula,
                          start = NULL,
                          dispersion = NULL,
                          family = stop("Family argument is needed"),
                          Fit.Method = "IRLS",
                          CI = "analytic",
                          prior.weights = 1,
                          ...) {
  if (is.character(family)){
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)){
    family <- family()
  }
  weights = 1
  Overdispersion <- dispersion
  Model_frame <- stats::model.frame(formula, data,  ...)
  Observed <- Model_frame[, 1]
  Variables <- stats::model.matrix(Model_frame)


  if (!family$valideta(start) && !is.null(start)){
    stop("Invalid start parameter")
  }

  log_like <- family$make_minusloglike(y = Observed, X = Variables, weight = prior.weights)
  grad <- family$make_gradient(y = Observed, X = Variables, weight = prior.weights)
  Hess <- family$make_hessian(y= Observed, X = Variables, weight = prior.weights)

  if (is.null(start)) {
    start <- stats::glm(formula = formula,
                 family = stats::poisson(),
                 data = data)$coefficients
  }

  if (family$family == "ZTP" || family$family == "ZTOIP"){
    df.reduced <- length(Observed) - dim(Variables)[2]
    if (Fit.Method == "IRLS") {
      FITT <- IRLS(Dependent = Observed,
                   Covariates = as.matrix(Variables),
                   eps = 1e-10,
                   family = family,
                   weights = prior.weights,
                   start = start)
      iter <- FITT$iter
      weights <- FITT$Weights
      coefficients <- FITT$Coefficients
    } else if (Fit.Method == "MaxLikelihood") {
      FITT <- stats::optim(fn = log_like,
                           par = start,
                           gr = grad)
      iter <- FITT$counts[1]
      coefficients <- FITT$par
    } else {
      stop("Method not implemented")
    }

    coefficients <- as.vector(coefficients)
    names(coefficients) <- colnames(Variables)
    Eta <- as.matrix(Variables) %*% coefficients
    Fitted <- family$linkinv(Eta)
    hess <- Hess(beta = coefficients)
  } else if (family$family == "ZTNB") {
    df.reduced <- length(Observed) - dim(Variables)[2]
    if (Fit.Method == "IRLS") {
      # stop("Method not implemented/n")
      # Overdispersion <- 0.5045529
      # TODO this is bad find a better estimator that doesn't depend on mean
      # this is method of moments estimate for non truncated neg-binomial
      # for truncated method of momnets lambert W function is needed
      # good luck future me you'll need it
      # Maxlikelihood works very well tho also if you already know overdispersion this works well
      if (is.null(Overdispersion)) {
        Overdispersion <- abs(mean(Observed ** 2) - mean(Observed)) / (mean(Observed) ** 2)
      }
      FITT <- IRLS(Dependent = Observed,
                   Covariates = as.matrix(Variables),
                   eps = 1e-10,
                   disp = Overdispersion,
                   family = Zero_Truncated_Negative_Binomial(),
                   weights = prior.weights,
                   start = start)

      iter <- FITT$iter
      weights <- FITT$Weights
      beta <- FITT$Coefficients

      coefficients <- c(Overdispersion, beta)
      Eta <- as.matrix(Variables) %*% beta
      Fitted <- family$linkinv(Eta)
    } else if (Fit.Method == "MaxLikelihood") {
      start <- c(.1, start)
      FITT <- stats::optim(fn = log_like,
                           par = start,
                           gr = grad)
      Beta <- FITT$par[-1]
      iter <- FITT$counts[1]
      Eta <- as.matrix(Variables) %*% Beta
      Fitted <- family$linkinv(Eta)

      C1 <- FITT$par[1]
      Overdispersion <- C1
      coefficients <- c(C1, Beta)
    } else {
      stop("Method not implemented/n")
    }

    names(coefficients) <- c("(Overdispersion)", colnames(Variables))
    hess <- Hess(ARG = coefficients)
  }

  Eta <- as.numeric(Eta)
  names(Eta) <- rownames(Variables)

  Std_err <- sqrt(diag(solve(-hess)))
  Tval <- coefficients / Std_err

  null.deviance <- as.numeric(NULL)
  LOG <- -log_like(coefficients)
  ResRes <- weights * (Observed - Fitted) #/ sqrt(family$variance(Eta))
  aic <- 2 * length(coefficients) - 2 * LOG
  deviance <- as.numeric(NULL)
  # VGAM::vglm używa przybliżenia rozkładem normalnym standaryzowanym można ewentualnie zamienić
  Pvals <- (stats::pt(q = abs(Tval), df = df.reduced, lower.tail = FALSE) +
            stats::pt(q = -abs(Tval), df = df.reduced, lower.tail = TRUE))

  Pop <- PopulationEstimate(y = Observed,
                            X = as.data.frame(Variables),
                            Grad = grad,
                            Hess = Hess,
                            Method = CI,
                            weights = prior.weights,
                            Parameter = Fitted,
                            family = family,
                            Overdispersion = Overdispersion,
                            beta = coefficients)


  Result <- list(y = Observed,
                 formula = formula,
                 call = match.call(),
                 coefficients = coefficients,
                 Standard_Errors = Std_err,
                 Tvalues = Tval,
                 Pvalues = Pvals,
                 null.deviance = null.deviance,
                 family = family,
                 aic = aic,
                 deviance = deviance,
                 prior.weights = prior.weights,
                 weights = weights,
                 residuals = ResRes,
                 LogL = LOG,
                 iter = iter,
                 df.residual = df.reduced,
                 df.null = length(Observed) - 1,
                 fitted.values = Fitted,
                 Population_Size = Pop,
                 model = Model_frame,
                 linear.predictors = Eta,
                 qr = qr(Variables))
  class(Result) <- c("glm", "lm")
  Result
}
