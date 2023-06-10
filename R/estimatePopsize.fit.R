#' \loadmathjax
#' 
#' @title Regression fitting in single source capture-recapture models
#' @author Piotr Chlebicki, Maciej Beresewicz
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
#' function will be used to fit regression with analytically computed gradient and 
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
#' 
#' 1. Initialise with:
#' \itemize{
#' \item \code{converged <- FALSE}
#' \item \code{iter <- 1}
#' \item \mjseqn{\boldsymbol{\beta}}\code{ <- start}
#' \item \mjseqn{\boldsymbol{W}}\code{ <- prior}
#' \item \mjseqn{\ell}\code{ <- }\mjseqn{\ell(\boldsymbol{\beta})}
#' }
#' 2. If \code{converged} or \code{iter > Maxiter} move to step 7.
#' 3. Store values from previous algorithm step:
#' \itemize{
#' \item \mjseqn{\boldsymbol{W}_{-}}\code{ <- }
#' \mjseqn{\boldsymbol{W}}
#' \item \mjseqn{\ell_{-}}\code{ <- }
#' \mjseqn{\ell}
#' \item \mjseqn{\boldsymbol{\beta}_{-}}\code{ <- }
#' \mjseqn{\boldsymbol{\beta}}
#' } and assign values at current step:
#' \itemize{
#' \item \mjseqn{\boldsymbol{\eta}}\code{ <- }
#' \mjseqn{\boldsymbol{X}_{vlm}\boldsymbol{\beta}}
#' \item \mjseqn{Z_{i}}\code{ <- }
#' \mjseqn{
#' \eta_{i}+\frac{\partial\ell_{i}}{\partial\eta_{i}}
#' \mathbb{E}\left(\frac{\partial^{2}\ell_{i}}{
#' \partial\eta_{i}^{T}\partial\eta_{i}}\right)^{-1}}
#' \item \mjseqn{\boldsymbol{W}_{ij}}\code{ <- }
#' \mjseqn{\mathbb{E}\left(\frac{\partial^{2}\ell}{
#' \partial\boldsymbol{\eta}_{j}^{T}\partial\boldsymbol{\eta}_{i}}\right)}
#' }
#' where \mjseqn{\ell_{i}} is the ith component of log likelihood 
#' function, \mjseqn{\eta_{i}} is the vector of linear predictors 
#' associated with ith row and \mjseqn{\mathbb{E}\left(\frac{\partial^{2}\ell_{i}}{
#' \partial\eta_{i}^{T}\partial\eta_{i}}\right)}
#' corresponds to weights associated with ith row and \mjseqn{\boldsymbol{W}} 
#' is a block matrix, made of diagonal matrixes 
#' \mjseqn{\mathbb{E}\left(\frac{\partial^{2}\ell}{
#' \partial\boldsymbol{\eta}_{j}^{T}\partial\boldsymbol{\eta}_{i}}\right)}
#' 4. Regress \mjseqn{\boldsymbol{Z}} on 
#' \mjseqn{\boldsymbol{X}_{vlm}} to obtain 
#' \mjseqn{\boldsymbol{\beta}} as:
#' \mjsdeqn{\boldsymbol{\beta}=
#' \left(\boldsymbol{X}_{vlm}^{T}\boldsymbol{W}\boldsymbol{X}_{vlm}\right)^{-1}
#' \boldsymbol{X}_{vlm}^{T}\boldsymbol{W}\boldsymbol{Z}}
#' 5. Assign:
#' \itemize{
#' \item\code{converged <- }\mjseqn{
#' \ell(\boldsymbol{\beta})-\ell_{-} < \varepsilon\cdot\ell_{-}}
#' or
#' \mjseqn{
#' ||\boldsymbol{\beta}-\boldsymbol{\beta}_{-}||_{\infty} < \varepsilon}
#' \item\code{iter <- iter + 1}
#' }
#' where \mjseqn{\varepsilon} is the relative tolerance level, 
#' by default \code{1e-8}.
#' 6. Return to step 2.
#' 7. Return \mjseqn{\boldsymbol{\beta}, \boldsymbol{W}}, \code{iter}.
#' 
#' In this package we use different conventions for \mjseqn{\boldsymbol{X}_{vlm}}
#' matrix hence slight differences are present in algorithm description but 
#' results are identical.
#' 
#' @references Yee, T. W. (2015). Vector Generalized Linear and Additive Models: 
#' With an Implementation in R. New York, USA: Springer. ISBN 978-1-4939-2817-0.
#' 
#' @examples 
#' \donttest{
#' # Get data
#' summary(farmsubmission)
#' 
#' # construct vglm model matrix
#' X <- matrix(data = 0, nrow = 2 * NROW(farmsubmission), ncol = 7)
#' X[1:NROW(farmsubmission), 1:4] <- model.matrix(
#'   ~ 1 + log_size + log_distance + C_TYPE, 
#'   farmsubmission
#' )
#' 
#' X[-(1:NROW(farmsubmission)), 5:7] <- X[1:NROW(farmsubmission), c(1, 3, 4)]
#' 
#' # this atrrtibute tells the function which elements of the design matrix 
#' # correspond to which linear predictor 
#' attr(X, "hwm") <- c(4, 3)
#' 
#' # get starting points
#' start <- glm.fit(
#'   y = farmsubmission$TOTAL_SUB, 
#'   x = X[1:NROW(farmsubmission), 1:4], 
#'   family = poisson()
#' )$coefficients
#'
#' start <- c(start, 0, 0, 0)
#' 
#' # call function
#' res <- estimatePopsize.fit(
#'   y = farmsubmission$TOTAL_SUB, 
#'   X = X, 
#'   method = "IRLS", 
#'   priorWeights = 1, 
#'   family = ztoigeom(), 
#'   control = controlMethod(verbose = 5), 
#'   start = start
#' )
#' 
#' # extract results
#' 
#' # regression coefficient vector
#' res$beta
#' 
#' # check likelihood
#' ll <- ztnegbin()$makeMinusLogLike(y = farmsubmission$TOTAL_SUB, X = X)
#' 
#' -ll(res$beta)
#' 
#' # number of iterations
#' res$iter
#' 
#' # working weights
#' head(res$weights)
#' 
#' # Compare with optim call
#' 
#' res2 <- estimatePopsize.fit(
#'   y = farmsubmission$TOTAL_SUB, 
#'   X = X, 
#'   method = "optim", 
#'   priorWeights = 1, 
#'   family = ztoigeom(), 
#'   start = start, 
#'   control = controlMethod(verbose = 1)
#' )
#' # extract results
#' 
#' # regression coefficient vector
#' res2$beta
#' 
#' 
#' # check likelihood
#' -ll(res2$beta)
#' 
#' # number of calls to log lik function
#' # since optim does not return the number of
#' # iterations
#' res2$iter
#' 
#' # optim does not calculated working weights
#' head(res2$weights)
#' }
#' @return List with regression parameters, working weights 
#' (if IRLS fitting method) was chosen and number of iterations taken.
#' @seealso [stats::glm()] [estimatePopsize()] [controlMethod()] [stats::optim()] 
# #' @importFrom maxLik maxLik
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
  logg <- NULL
  
  if (method == "IRLS") {
    
    FITT <- singleRcaptureinternalIRLSmultipar(
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
      crit = control$criterion,
      saveLog = control$saveIRLSlogs,
      printOften = control$printEveryN
    )
    
    iter <- FITT$iter
    weights <- FITT$weights
    beta <- FITT$coefficients
    logg <- FITT$logg
  } else if (method == "optim") {
    logLike <- family$makeMinusLogLike(
      y = y, X = X, 
      weight = priorWeights
    )
    grad <- family$makeMinusLogLike(
      y = y, X = X, 
      weight = priorWeights, 
      deriv = 1
    )
    
    weights <- priorWeights
    methodopt <- control$optimMethod
    
    if (!isFALSE(control$optimPass)) {
      ctrl <- control$optimPass
    } else {
      if (methodopt == "L-BFGS-B") {
        ctrl <- list(
          factr = control$epsilon,
          maxit = control$maxiter,
          trace = control$verbose
        )
      } else {
        ctrl <- list(
          reltol = control$epsilon,
          maxit = control$maxiter,
          trace = control$verbose
        )
      }
    }
    
    
    giveError <- tryCatch(
      expr = {
        FITT <- stats::optim(
          fn = logLike,
          par = start,
          gr = function(x) -grad(x),
          method = methodopt,
          # hessian = TRUE,
          control = ctrl
        )
        list(0)
      },
      error = function (e) {
        list(1, e)
      },
      warning = function (w) {
        list(2, w)
      }
    )

    if (giveError[[1]] == 1) {
      stop(paste0(
        "stats::optim fitting returned with the following error: ",
        giveError[[2]]$message,
        " consider fitting with 'IRLS' or changing controlMethod argument."
      ))
    }
    
    if (!isTRUE(control$silent) && (giveError[[1]] == 2)) {
      warning("stats::optim returned with the following warning: ",
              giveError[[2]]$message)
    }
    
    if (FITT$convergence != 0 && !isTRUE(control$silent)) {
      warning("Convergence not obtained in ", control$maxiter, 
              " iterations of optim fitting algorithm", sep = "")
      
      switch (as.character(FITT$convergence),
        "1"  = warning("Convergence not obtained in ", control$maxiter, 
                       " iterations of optim fitting algorithm", sep = ""),
        "10" = warning("stats::optim returned with code 10 degeneracy of the Nelder-Mead simplex."),
        "51" =    stop("stats::optim returned with code 51 with message: ", FITT$message),
        "52" = warning("stats::optim returned with code 52 with message: ", FITT$message)
      )
    }
    
    iter <- FITT$counts
    beta <- FITT$par
  } else if (method == "maxLik") {
    # weights <- priorWeights
    # ll <- family$makeMinusLogLike(y = y, X = X, weight = priorWeights)
    # gr <- family$makeMinusLogLike(y = y, X = X, weight = priorWeights, deriv = 1)
    # FITT <- maxLik::maxLik(
    #   logLik = function(x) -ll(x),
    #   grad = function(x) matrix(gr(x), nrow = 1),
    #   hess = family$makeMinusLogLike(y = y, X = X, weight = priorWeights, deriv = 2),
    #   method = "NM",
    #   start = start,
    #   iterlim = 10000
    # )
    # print(FITT)
    stop("Currently in development")
  } else {
    stop("Method implemented.")
  }
  if (is.null(FITT)) {
    stop("fitting error try another model
          (negative binomial models are highly volitile)")
  }

  beta <- as.vector(beta)

  list(
    beta = beta,
    weights = weights,
    iter = iter,
    logg = logg
  )
}
