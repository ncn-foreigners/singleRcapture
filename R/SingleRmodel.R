#' SingleRmodel
#'
#' @param data Data frame for the regression
#' @param Dependent Name of column of dependent variables in data
#' @param Explanatory Vector of names of columns of explanatory variables (covariates) in data
#' @param Fit.Method Method for fitting values currently supported IRLS and MaxLikelihood
#' @param Family Family of distributions used in regression
#'
#' @return Returns a S3 object with:\cr fitted values\cr standard errors\cr akane information criterion \cr
#' Tvalues and Pvalues for testing the hypothesis of no effect of coefficient on expected value \cr
#' Names of used variables\cr response residuals\cr estimate of population size pearson residuals and parameter\cr
#' Estimation of true population size
#' @export
SingleRmodel <- function(data,Dependent,
                                    Explanatory,Family = stop('Family is needed'),Fit.Method = 'IRLS'){
  Observed <- data[,Dependent]
  Variables <- data[,Explanatory]
  Variables <- cbind(Intercept = 1,Variables)

  log_like <- Family()$make_minusloglike(y = Observed,X = Variables)
  grad <- Family()$make_gradient(y = Observed,X = Variables)
  Hess <- Family()$make_hessian(X = Variables)

  if(Fit.Method == 'IRLS'){
    Predictors <- IRLS(Observed,as.matrix(Variables),eps = 1e-20,
                                      Family = Family,start = rep(.1,length(Variables)))$Coefficients
  }

  else if(Fit.Method == 'MaxLikelihood'){
    Predictors <- stats::optim(fn = log_like, par = rep(.1,length(Variables)),gr = grad,method = 'L-BFGS-B')
  }
  else{
    stop('Method not implemented/n')
  }
  names(Predictors) <- colnames(Variables)

  hess <- Hess(Predictors)
  Std_err <- sqrt(diag(solve(-hess)))
  Tval <- Predictors / Std_err

  Parameter <- exp(as.matrix(Variables) %*% Predictors)
  Expected <- Family()$mu(Parameter)

  ResRes <- Observed-Expected
  PearsonRes <- (Observed-Expected) / stats::sd(Expected)
  AIC <- 2 * (1+log_like(Predictors))
  df.reduced <- length(Observed) - length(Variables)
  # VGLM używa przybliżenia rozkładem normalnym można ewentualnie zamienić
  Pvals <- stats::pt(q = abs(Tval),df = df.reduced,lower.tail = F)+stats::pt(q = -abs(Tval),df = df.reduced,lower.tail = T)

  Pop <- PoissonPopEst(y = Observed, X = as.data.frame(Variables) , Grad = grad,
                 Hess = Hess,lambda = Parameter, beta = Predictors)

  # Deviance residuals
  Result <- list(Coefficients = Predictors,StandardErrors = Std_err,Tvalues = Tval,Pvalues = Pvals,
                 ResponseResiduals = ResRes,ExplanatoryVariables=Explanatory,DependentVariable = Dependent,AIC = AIC,
                 LogL = -log_like(Predictors),df.reduced = df.reduced,Parameter = Parameter,
                 PearsonResiduals = PearsonRes,PopulationSize = Pop)
  # Figure out a better name for the class
  class(Result) <- 'Model'
  return(Result)
}
