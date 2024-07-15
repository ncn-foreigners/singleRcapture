#' @title Generating data in singleRcapture
#' 
#' @description
#' An S3 method for \code{stats::simulate} to handle \code{singleRStaticCountData} and 
#' \code{singleRfamily} classes.
#'
#' @param object an object representing a fitted model.
#' @param nsim a numeric scalar specifying:
#' \itemize{
#'    \item number of response vectors to simulate in \code{simulate.singleRStaticCountData}, defaults to \code{1L}.
#'    \item number of units to draw in \code{simulate.singleRfamily}, defaults to \code{NROW(eta)}.
#' }
#' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).
#' @param truncated logical value indicating whether to sample from truncated or
#' full distribution.
#' @param eta a matrix of linear predictors
#' @param ... additional optional arguments.
#' @return a \code{data.frame} with \code{n} rows and \code{nsim} columns.
#' @seealso [stats::simulate()] [singleRcapture::estimatePopsize()]
#' @examples
#' N <- 10000
#' ###gender <- rbinom(N, 1, 0.2)
#' gender <- rep(0:1, c(8042, 1958))
#' eta <- -1 + 0.5*gender
#' counts <- simulate(ztpoisson(), eta = cbind(eta), seed = 1)
#' df <- data.frame(gender, eta, counts)
#' df2 <- subset(df, counts > 0)
#' ### check coverage with summary
#' mod1 <-  estimatePopsize(
#'   formula       = counts ~ 1 + gender, 
#'   data          = df2, 
#'   model         = ztpoisson, 
#'   controlMethod = list(silent = TRUE)
#' )
#' mod1_sims <- simulate(mod1, nsim=10, seed = 1)
#' colMeans(mod1_sims)
#' mean(df2$counts)
#' @author Maciej Beręsewicz, Piotr Chlebicki
#' @importFrom stats simulate
#' @method simulate singleRStaticCountData
#' @exportS3Method
#' @name simulate
#' @export
simulate.singleRStaticCountData <- function(object, nsim = 1, seed = NULL, ...) {
  n <- nobs(object)
  eta <- object$linearPredictors
  
  # Replicate each row in eta priorWeights number of times
  if (isTRUE(object$control$controlModel$weightsAsCounts)) {
    eta <- matrix(
      unlist(sapply(1:NROW(eta), function(x) {
        rep(eta[x,], object$priorWeights[x])
      })),
      ncol = NCOL(eta),
      byrow = TRUE,
      dimnames = list(
        1:n, colnames(eta)
      )
    )
  }
  
  val <- simulate(
    object    = family(object), 
    seed      = seed, 
    nsim      = n * nsim, 
    eta       = eta, 
    truncated = TRUE
  )

  dim(val) <- c(n, nsim)
  val <- as.data.frame(val)
  names(val) <- paste0("sim_", seq_len(nsim))
  val
}

#' @rdname simulate
#' @importFrom stats simulate
#' @method simulate singleRfamily
#' @exportS3Method
#' @export
simulate.singleRfamily <- function(object, 
                                   nsim,
                                   seed = NULL, 
                                   eta, 
                                   truncated = FALSE, 
                                   ...) {
  if (missing(nsim)) nsim <- NROW(eta)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  object$simulate(nsim, eta, lower = ifelse(truncated, 0, -1))
}
