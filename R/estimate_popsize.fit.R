#' @title Regression fitting in single source capture-recapture models
#'
#' @description \code{estimate_popsize.fit} does for \code{estimate_popsize} what
#' \code{glm.fit} does for \code{glm}. It is internally called in 
#' \code{estimate_popsize}. Since \code{estimate_popsize} does much more than
#' just regression fitting \code{estimate_popsize.fit} is much faster.
#'
#' @param y vector of dependent variables.
#' @param X model matrix, the vglm one.
#' @param family same as model in \code{estimate_popsize}.
#' @param control control parameters created in \code{control.model}.
#' @param method method of estimation same as in \code{estimate_popsize}.
#' @param prior.weights vector of prior weights its the same argument as weights
#' in \code{estimate_popsize}.
#' @param start initial value of regression parameters.
#' @param ... arguments to pass to other methods.
#' 
#' @details If \code{method} argument was set to \code{"mle"} the \code{stats::optim}
#' function will be used to fit regression with analyticly computed gradient and 
#' (minus) log likelihood functions as \code{gr} and \code{fn} arguments. 
#' Unfortunately \code{optim} does not allow for hessian to be specified.
#' More information about how to modify \code{optim} fitting is included in 
#' [control.method()].
#' 
#' 
#' If \code{method} argument was set to \code{"robust"} the iteratively reweighted 
#' least squares. The algorithm is well know in simple generalised linear models.
#' Thomas W. Yee later extended this algorithm to vector generalised linear models 
#' and in more general terms it can roughly be described as:
#' \loadmathjax
#' 1. Initialise with:
#' \itemize{
#' \item \code{converged <- FALSE}
#' \item \code{iter <- 1}
#' \item \mjeqn{\boldsymbol{\beta}}{beta}\code{ <- start}
#' \item \mjeqn{\boldsymbol{W}}{W}\code{ <- prior}
#' \item \mjeqn{\ell}{l}\code{ <- }\mjeqn{\ell(\boldsymbol{\beta})}{l(beta)}
#' }
#' 2. If \code{converged} or \code{iter > Maxiter} move to step 7
#' 3. Store values from previous algorithm step:
#' \itemize{
#' \item \mjeqn{\boldsymbol{W}_{-}}{W_-}\code{ <- }\mjeqn{\boldsymbol{W}}{W}
#' \item \mjeqn{\ell_{-}}{l_-}\code{ <- }\mjeqn{\ell}{l}
#' \item \mjeqn{\boldsymbol{\beta}_{-}}{beta_-}\code{ <- }\mjeqn{\boldsymbol{\beta}}{beta}
#' } and assign values at current step:
#' \itemize{
#' \item \mjeqn{\boldsymbol{\eta}}{eta}\code{ <- }\mjeqn{\boldsymbol{X}_{vlm}\boldsymbol{\beta}}{X_vlm * beta}
#' \item \mjeqn{Z_{i}}{Z_i}\code{ <- }\mjeqn{\eta_{i}+\frac{\partial\ell_{i}}{\partial\eta_{i}}\mathbb{E}\left(\frac{\partial^{2}\ell_{i}}{\partial\eta_{i}^{T}\partial\eta_{i}}\right)^{-1}}{eta_i + (dl/deta_i)(E(d^2l/deta_i^Tdeta_i))^-1}
#' \item \mjeqn{\boldsymbol{W}_{ij}}{W_ij}\code{ <- }\mjeqn{\mathbb{E}\left(\frac{\partial^{2}\ell}{\partial\boldsymbol{\eta}_{j}^{T}\partial\boldsymbol{\eta}_{i}}\right)}{E(d^2l/deta_j^Tdeta_i)}
#' }
#' where \mjeqn{\ell_{i}}{l_i} is the ith component of log likelihood function, 
#' \mjeqn{\eta_{i}}{eta_i} is the vector of linear predictors associated with 
#' ith row and  \mjeqn{\mathbb{E}\left(\frac{\partial^{2}\ell_{i}}{\partial\eta_{i}^{T}\partial\eta_{i}}\right)}{E(d^2l_i/deta_i^Tdeta_i)} 
#' corresponds to weights associated with ith row and \mjeqn{\boldsymbol{W}}{W} 
#' is a block matrix.
#' 4. Regress \mjeqn{\boldsymbol{Z}}{Z} on \mjeqn{\boldsymbol{X}_{vlm}}{X_vlm} 
#' to obtain \mjeqn{\boldsymbol{\beta}}{beta} as:
#' \mjdeqn{\boldsymbol{\beta}=\left(\boldsymbol{X}_{vlm}^{T}\boldsymbol{W}\boldsymbol{X}_{vlm}\right)^{-1}\boldsymbol{X}_{vlm}^{T}\boldsymbol{W}\boldsymbol{Z}}{(X_vlm^T W X_vlm)^-1(X_vlm^T W Z)}
#' 5. Assign:
#' \itemize{
#' \item\code{converged <- }\mjeqn{\ell(\boldsymbol{\beta})-\ell_{-} < \ell_{-}\cdot\varepsilon}{l(beta)-l_-<l_- * epsilon} \code{or} 
#' \mjeqn{||\boldsymbol{\beta}-\boldsymbol{\beta}_{-}||_{\infty} < \varepsilon}{||beta - beta_-||_inf< epsilon}
#' \item\code{iter <- iter + 1}
#' }
#' where \mjeqn{\varepsilon}{epsilon} is the relative tolerance level, by default \code{1e-8}.
#' 6. Return to step 2
#' 7. Return \mjeqn{\boldsymbol{\beta}, \boldsymbol{W}}{beta, W}, \code{iter}
#' 
#' @references Yee, T. W. (2015). Vector Generalized Linear and Additive Models: 
#' With an Implementation in R. New York, USA: Springer. ISBN 978-1-4939-2817-0.
#' 
#' @return List with regression parameters, working weights 
#' (if IRLS fitting method) was chosen and number of iterations taken.
#' @seealso [stats::glm()] [estimate_popsize()] [control.method()] [stats::optim()]
#' @export
estimate_popsize.fit <- function(y, X,
                                 family,
                                 control,
                                 method,
                                 prior.weights,
                                 start,
                                 ...) {
  hwm <- attr(X, "hwm")
  tbgname <- colnames(X)
  X <- as.matrix(X)

  if (grepl("hurdle", family$family)) {
    FITT <- singleRcaptureinternalIRLS(dependent = y,
                                       covariates = X,
                                       eps = control$epsilon,
                                       family = singleRcapture::chao(),
                                       maxiter = control$maxiter,
                                       weights = prior.weights,
                                       start = start,
                                       silent = control$silent,
                                       trace = control$verbose,
                                       stepsize = control$stepsize)
    
    gamma <- FITT$coefficients
    
    FITT <- singleRcaptureinternalIRLS(dependent = y,
                                       covariates = X,
                                       eps = control$epsilon,
                                       family = chao(),
                                       maxiter = control$maxiter,
                                       weights = prior.weights,
                                       start = start,
                                       silent = control$silent,
                                       trace = control$verbose,
                                       stepsize = control$stepsize)
    
    weights <- FITT$weights
    beta <- FITT$coefficients
  } else {
    if (method == "robust") {
      
      if (family$parNum == 1) {
        FittingFunction <- singleRcaptureinternalIRLS
      } else {
        FittingFunction <- singleRcaptureinternalIRLSmultipar
      }
      
      FITT <- FittingFunction(
        dependent = y,
        covariates = X,
        eps = control$epsilon,
        family = family,
        maxiter = control$maxiter,
        weights = prior.weights,
        start = start,
        silent = control$silent,
        trace = control$verbose,
        stepsize = control$stepsize,
        hwm = hwm,
        momentumFactor = control$momentumFactor,
        momentumActivation = control$momentumActivation,
      )
      
      iter <- FITT$iter
      weights <- FITT$weights
      beta <- FITT$coefficients
    } else if (method == "mle") {
      logLike <- family$makeMinusLogLike(y = y, X = X, weight = prior.weights)
      grad <- family$makeGradient(y = y, X = X, weight = prior.weights, lambdaPredNumber = length(start) - 1)
      
      weights <- prior.weights
      methodopt <- control$mleMethod
      
      if (!isFALSE(control$optimPass)) {
        ctrl <- control$optimPass
      } else {
        if (methodopt == "L-BFGS-B") {
          ctrl <- list(
            factr = control$epsilon,
            maxit = control$maxiter,
            trace = control$verbose,
            trace = if (is.numeric(control$trace)) control$trace else 0
          )
        } else {
          ctrl <- list(
            reltol = control$epsilon,
            maxit = control$maxiter,
            trace = control$verbose
          )
        }
      }

      FITT <- stats::optim(fn = logLike,
                           par = start,
                           gr = function(x) -grad(x),
                           method = methodopt,
                           # hessian = TRUE,
                           control = ctrl)
      beta <- FITT$par
      # print(FITT)
      # print(rootSolve::gradient(f = grad, x = beta))
      # Commented lines of code are used in verification of computed analytic hessian
      iter <- FITT$counts
      if (FITT$convergence == 1 && !control$silent) {
        warning("Convergence not obtained in ", control$maxiter, " iterations of mle fitting algorithm", sep = "")
      }
    } else {
      stop("Method not implemented")
    }
  }

  if (is.null(FITT)) {
    stop("fitting error try another model
          (negative binomial models are highly volitile)")
  }

  beta <- as.vector(beta)

  list(
    beta = beta,
    weights = weights,
    iter = iter
  )
}
