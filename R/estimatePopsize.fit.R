#' @title Regression fitting in single source capture-recapture models
#' @author Piotr Chlebicki, Maciej Beręsewicz
#'
#' @description \code{estimatePopsize.fit} does for \code{estimatePopsize} what
#' \code{glm.fit} does for \code{glm}. It is internally called in 
#' \code{estimatePopsize}. Since \code{estimatePopsize} does much more than
#' just regression fitting \code{estimatePopsize.fit} is much faster.
#'
#' @param y vector of dependent variables.
#' @param X model matrix, the vglm one.
#' @param family same as model in \code{estimatePopsize}.
#' @param control control parameters created in \code{controlModel}.
#' @param method method of estimation same as in \code{estimatePopsize}.
#' @param priorWeights vector of prior weights its the same argument as weights
#' in \code{estimatePopsize}.
#' @param start initial value of regression parameters.
#' @param ... arguments to pass to other methods.
#' 
#' @details If \code{method} argument was set to \code{"optim"} the \code{stats::optim}
#' function will be used to fit regression with analyticly computed gradient and 
#' (minus) log likelihood functions as \code{gr} and \code{fn} arguments. 
#' Unfortunately \code{optim} does not allow for hessian to be specified.
#' More information about how to modify \code{optim} fitting is included in 
#' [controlMethod()].
#' 
#' 
#' If \code{method} argument was set to \code{"IRLS"} the iteratively reweighted 
#' least squares. The algorithm is well know in generalised linear models.
#' Thomas W. Yee later extended this algorithm to vector generalised linear models 
#' and in more general terms it can roughly be described as 
#' (this is Yee's description after changing some conventions):
#' \loadmathjax
#' 1. Initialise with:
#' \itemize{
#' \item \code{converged <- FALSE}
#' \item \code{iter <- 1}
#' \item \mjeqn{\boldsymbol{\beta}}{beta}\code{ <- start}
#' \item \mjeqn{\boldsymbol{W}}{W}\code{ <- prior}
#' \item \mjeqn{\ell}{l}\code{ <- }\mjeqn{\ell(\boldsymbol{\beta})}{l(beta)}
#' }
#' 2. If \code{converged} or \code{iter > Maxiter} move to step 7.
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
#' is a block matrix, made of diagonal matrixes \mjeqn{\mathbb{E}\left(\frac{\partial^{2}\ell}{\partial\boldsymbol{\eta}_{j}^{T}\partial\boldsymbol{\eta}_{i}}\right)}{E(d^2l/deta_j^Tdeta_i)}.
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
#' 6. Return to step 2.
#' 7. Return \mjeqn{\boldsymbol{\beta}, \boldsymbol{W}}{beta, W}, \code{iter}.
#' 
#' In this package we use different conventions for \mjeqn{\boldsymbol{X}_vlm}{X_vlm}
#' matrix hence slight differences are present in algorithm description but
#' results are identical.
#' 
#' @references Yee, T. W. (2015). Vector Generalized Linear and Additive Models: 
#' With an Implementation in R. New York, USA: Springer. ISBN 978-1-4939-2817-0.
#' 
#' @examples 
#' 
#' # Set seed for reproducibility
#' set.seed(123)
#' 
#' # Set population size and 'true' regression coefficients
#' N <- 7562
#' beta1 <- c(1, .2, -.95, -3)
#' beta2 <- c(2, 0, -.01, .3)
#' 
#' # Generate covariates
#' x1 <- rpois(n = N, lambda = 1.85)
#' x2 <- rnorm(n = N)
#' x3 <- rbinom(n = N, size = 1, prob = runif(n = 1))
#' 
#' # Compute linear predictors
#' eta1 <- cbind(1, x1, x2, x3) %*% beta1
#' eta2 <- cbind(1, x1, x2, x3) %*% beta2
#' 
#' # generate dependent variable
#' y <- ztoipoisson()$simulate(n = N, eta = cbind(eta1, eta2))
#' # optionally visualise destribution on y with barplot
#' # barplot(table(y))
#' # truncate data
#' x1 <- x1[y>0]
#' x2 <- x2[y>0]
#' x3 <- x3[y>0]
#' y <- y[y>0]
#' 
#' # construct vglm model matrix
#' X <- matrix(data = 0, nrow = 2 * length(y), ncol = 8)
#' X[1:length(y), 1:4] <- model.matrix(~ 1 + x1 + x2 + x3)
#' X[-(1:length(y)), 5:8] <- X[1:length(y), 1:4]
#' # this arrtibute tells the function which elements of the design matrix 
#' # correspond to which linear predictor 
#' attr(X, "hwm") <- c(4, 4)
#' 
#' # finally call the function
#' beta <- c(beta1, beta2)
#' 
#' res <- estimatePopsize.fit(
#'   y = y, X = X, 
#'   method = "IRLS", 
#'   priorWeights = 1, 
#'   family = ztoipoisson(), 
#'   control = controlMethod(verbose = 5), 
#'   start = beta + rnorm(n = length(beta), sd = .1)
#' )
#' 
#' # extract results
#' 
#' # regression coefficient vector
#' res$beta
#' 
#' # compare with actual true coeffictient vector
#' beta
#' 
#' # number of iterations
#' res$iter
#' 
#' # working weights
#' head(res$weights)
#' 
#' # call with another method
#' 
#' res <- estimatePopsize.fit(
#'   y = y, X = X, 
#'   method = "optim", 
#'   priorWeights = 1, 
#'   family = ztoipoisson(), 
#'   start = beta + rnorm(n = length(beta), sd = .1), 
#'   control = controlMethod()
#' )
#' # extract results
#' 
#' # regression coefficient vector
#' res$beta
#' 
#' # compare with actual true coeffictient vector
#' beta
#' 
#' # number of iterations
#' res$iter
#' 
#' # optim does not calculated working weights
#' head(res$weights)
#' 
#' @return List with regression parameters, working weights 
#' (if IRLS fitting method) was chosen and number of iterations taken.
#' @seealso [stats::glm()] [estimatePopsize()] [controlMethod()] [stats::optim()]
#' @export
estimatePopsize.fit <- function(y, X,
                                family,
                                control,
                                method,
                                priorWeights,
                                start,
                                ...) {
  hwm <- attr(X, "hwm")
  tbgname <- colnames(X)
  X <- as.matrix(X)

  if (FALSE) {
    # maybe implement
  } else {
    if (method == "IRLS") {
      
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
        weights = priorWeights,
        start = start,
        silent = control$silent,
        trace = control$verbose,
        stepsize = control$stepsize,
        hwm = hwm,
        momentumFactor = control$momentumFactor,
        momentumActivation = control$momentumActivation,
        check = control$checkDiagWeights,
        epsWeights = control$weightsEpsilon,
        crit = control$criterion
      )
      
      iter <- FITT$iter
      weights <- FITT$weights
      beta <- FITT$coefficients
    } else if (method == "optim") {
      logLike <- family$makeMinusLogLike(y = y, X = X, weight = priorWeights)
      grad <- family$makeMinusLogLike(y = y, X = X, weight = priorWeights, deriv = 1)
      
      weights <- priorWeights
      methodopt <- control$optimMethod
      
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
        warning("Convergence not obtained in ", control$maxiter, " iterations of optim fitting algorithm", sep = "")
      }
    } else if (method == "maxLik") {
      #TODO
      stop("Not yet implemented")
    } else {
      stop("Method implemented.")
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
