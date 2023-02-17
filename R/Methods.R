#' @title Summary statistics for model of singleR class.
#' 
#' @description A \code{summary} method for \code{singleR} class, works 
#' analogically to \code{summary.glm} but includes population size estimation 
#' results. If any additional statistics, such as confidence intervals for
#' coefficients or coefficient correlation, are specified they will be printed.
#' 
#' @param object Object of singleR class.
#' @param test Type of test for significance of parameters \code{"t"} for t-test 
#' and \code{"z"} for normal approximation of students t distribution, by 
#' default \code{"z"} is used if there are more than 30 degrees of freedom
#' and \code{"t"} is used in other cases.
#' @param resType Type of residuals to summarise any value that is allowed in
#' \code{residuals.signleR} except for \code{"all"} is allowed. By default 
#' pearson residuals are used.
#' @param correlation Logical value indicating whether correlation matrix should
#' be computed from covariance matrix by default \code{FALSE}.
#' @param confint Logical value indicating whether confidence intervals for
#' regression parameters should be constructed. By default \code{FALSE}.
#' @param cov Covariance matrix corresponding to regression parameters. 
#' It is possible to give \code{cov} argument as a function of \code{object}.
#' If not specified it will be constructed using \code{vcov.singleR} method.
#' (i.e using Cramer-Rao lower bound)
#' @param popSizeEst A \code{popSizeEstResults} class object.
#' If not specified population size estimation results will be drawn from
#' \code{object}. If any post-hoc procedures, such as sandwich covariance matrix 
#' estimation or bias reduction, were taken it is possible to include them in 
#' population size estimation results by calling \code{redoPopEstimation}.
#' @param ... Additional optional arguments passed to the following functions:
#' \itemize{
#' \item \code{vcov.singleR} -- if no \code{cov} argument was provided.
#' \item \code{cov} -- if \code{cov} parameter specified at call was a function.
#' \item \code{confint.singleR} -- if \code{confint} parameter was set to \code{TRUE} at function call.
#' In particular it is possible to set confidence level in \code{...}.
#' }
#' 
#' @return An object of \code{summarysingleR} class containing:
#' \itemize{
#' \item \code{call} -- A call which created \code{object}.
#' \item \code{coefficients} -- A dataframe with estimated regression coefficients
#' and their summary statistics such as standard error wald test statistic and
#' p value for wald test.
#' \item \code{residuals} -- A vector of residuals of type specified at call.
#' \item \code{aic} -- Akaike's informsation criterion.
#' \item \code{bic} -- Bayesian (Schwarz's) information criterion.
#' \item \code{iter} -- Number of iterations taken in fitting regression.
#' \item \code{logL} -- Logarithm of likelihood function evaluated at coefficients.
#' \item \code{deviance} -- Residual deviance.
#' \item \code{populationSize} -- Object with population size estimation results.
#' \item \code{dfResidual} -- Residual degrees of freedom.
#' \item \code{sizeObserved} -- Size of observed population.
#' \item \code{correlation} -- Correlation matrix if \code{correlation} parameter was set to \code{TRUE}
#' \item \code{test} -- Type of statistical test performed.
#' \item \code{model} -- Family class object specified in call for \code{object}.
#' \item \code{skew} -- If bootstrap sample was saved contains estimate of skewness.
#' }
#' @method summary singleR
#' @importFrom stats pt
#' @importFrom stats coef
#' @importFrom stats sd
#' @seealso [redoPopEstimation()] [stats::summary.glm()]
#' @exportS3Method 
summary.singleR <- function(object, 
                            test = c("t", "z"), 
                            resType = "pearson", 
                            correlation = FALSE, 
                            confint = FALSE, 
                            cov, 
                            popSizeEst, 
                            ...) {
  if (resType == "all") {stop("Can't use 'resType = all' in summary.singleR method, if you wish to obtain all aviable types of residuals call residuals.singleR method directly.")}
  dfResidual <- object$dfResidual
  if (missing(test)) {if (dfResidual > 30) test <- "z" else test <- "t"}
  if (missing(cov)) {
    cov <- vcov.singleR(object, ...)
  } else if (is.function(cov)) {
    cov <- cov(object, ...)
  }
  pers <- residuals.singleR(object, type = resType)
  cf <- coef(object)
  se <- sqrt(diag(cov))
  wValues <- cf / se
  pValues <- switch (test,
  "t" = 2 *    stats::pt(q = -abs(wValues), df = dfResidual),
  "z" = 2 * stats::pnorm(q = abs(wValues), lower.tail = FALSE)
  )
  crr <- if (isFALSE(correlation)) {NULL} else {cov / outer(se, se)}
  if(isTRUE(correlation)) {rownames(crr) <- colnames(crr) <- names(cf)}
  cnfint <- if(isTRUE(confint)) {confint.singleR(object, ...)} else {NULL}
  if (is.numeric(object$populationSize$boot)) {
    n <- length(object$populationSize$boot)
    m <- sum((object$populationSize$boot - mean(object$populationSize$boot)) ^ 3) / n
    s <- sd(object$populationSize$boot)
    skew <- m / (s ^ 3)
  } else {
    skew <- NULL
  }
  
  ob <- data.frame(cf, se, wValues, pValues)
  
  colnames(ob) <- switch(test,
  "t" = c("Estimate", "Std. Error", "t value", "P(>|t|)"),
  "z" = c("Estimate", "Std. Error", "z value", "P(>|z|)"))
  if (isTRUE(confint)) {
    ob[, 4] <- cnfint[, 1]
    ob[, 5] <- cnfint[, 2]
    ob[, 6] <- pValues
    colnames(ob)[4:6] <- c(colnames(cnfint), colnames(ob)[4])
  }
  
  cf <- ob
  structure(
    list(
      call = object$call,
      coefficients = cf,
      residuals = pers,
      aic = AIC(object, ...),
      bic = BIC(object, ...),
      iter = object$iter,
      logL = object$logL,
      deviance = object$deviance,
      populationSize = if (missing(popSizeEst)) object$populationSize else popSizeEst,
      dfResidual = dfResidual,
      sizeObserved = object$sizeObserved,
      correlation = crr,
      test = test,
      model = object$model,
      skew = skew
    ),
    class = "summarysingleR"
  )
}
#' @title Regression deletion diagnostics for Population size estimation.
#'
#' @description \code{dfpopsize} does for population size estimation.
#'
#' @param model Model for which leave one out diagnostic of popsize will be done.
#' @param dfbeta If dfbeta was already obtained it is possible to pass them into 
#' function so that they need not be computed for the second time.
#' @param observedPop Logical. For \code{singleR} class object if set to 
#' \code{TRUE} indicates that 1 will be returned for units which do not
#' take part in population size estimation (e.g. 1's in zero one truncated
#' models or units with count => 3 for \code{zelterman} of basic \code{chao}
#' model) if set to \code{FALSE} (default) these units will not be included
#' in results.
#' @param ... Additional optional arguments passed to the following functions 
#' (for \code{singleR} class):
#' \itemize{
#' \item dfbetas -- if \code{dfbeta} parameter was not provided.
#' }
#'
#' @return A vector for which k'th element corresponds to the difference between
#' point estimate of population size estimation on full data set and 
#' point estimate of population size estimation after the removal of k'th
#' unit from the data set.
#' @examples 
#' \dontrun{
#' # For singleR class
#' # Get simple model
#' Model <- estimatePopsize(formula = capture ~ nation + age + gender, 
#' data = netherlandsimmigrant, 
#' model = ztpoisson, 
#' method = "IRLS")
#' # Get df beta
#' dfb <- dfbeta(Model)
#' # The results
#' dfpopsize(Model, dfbeta = dfb)
#' # It is also possible to not provide dfbeta then they will be
#' # computed manually
#' dfpopsize(Model)
#' }
#' @export
dfpopsize <- function(model, ...) {
  UseMethod("dfpopsize")
}
#' @title Updating population size estimation results.
#'
#' @description A function that applies all post-hoc procedures that were taken
#' (such as heteroscedastic consistent covariance matrix estimation or bias
#' reduction) to population size estimation and standard error estimation.
#' 
#' @param object Object for which update of population size estimation results will be done.
#' @param cov An updated covariance matrix estimate.
#' @param ... Additional optional arguments, currently not used in \code{singleR} class method.
#'
#' @return An object of class \code{popSizeEstResults} containing updated 
#' population size estimation results.
#' @examples
#' # Create simple model
#' Model <- estimatePopsize(formula = capture ~ nation + gender, 
#' data = netherlandsimmigrant, 
#' model = ztpoisson, 
#' method = "IRLS")
#' # Apply heteroscedasticity consistent covariance matrix estimation
#' require(sandwich)
#' cov <- vcovHC(Model, type = "HC3")
#' summary(Model, cov = cov,
#' popSizeEst = redoPopEstimation(Model, cov = cov))
#' # Compare to results with usual covariance matrix estimation
#' summary(Model)
#' @export
redoPopEstimation <- function(object, ...) {
  UseMethod("redoPopEstimation")
}
#' @title Extract population size estimation results.
#' 
#' @description An extractor function with \code{singleR} method for extracting
#' important information regarding pop size estimate.
#'
#' @param object Object with population size estimates.
#' @param ... Additional optional arguments, currently not used in \code{singleR} class method. 
#'
#' @return An object of class \code{popSizeEstResults} containing population size estimation results.
#' @export
popSizeEst <- function(object, ...) {
  UseMethod("popSizeEst")
}
#' @title Estimate size of sub populations.
#' 
#' @description A function that estimates sizes of specific sub populations 
#' based on a capture-recapture model for the whole population.
#'
#' @param object An object on which the population size estimates should be based
#' in \code{singleRcapture} package this is a fitter \code{singleR} class object.
#' @param stratas A specification of sub populations either by:
#' \itemize{
#' \item formula -- TODO.
#' \item logical vector with number of entries equal to number of rows in the dataset.
#' \item A (named) list where each element is a logical vector, names of the list
#' will be used to specify row names in returned object.
#' \item Vector of names of explanatory variables. For \code{singleR} method
#' for this function this specification of \code{stratas} parameter will
#' result in every level of explanatory variable having its own sub population
#' for each variable specified.
#' \item If no value was provided the \code{singleR} method for this function 
#' will itself create sub populations based on levels of factor variables
#' in \code{model.frame}.
#' }
#' @param cov For \code{singleR} method an estimate of variance-covariance matrix
#' for estimate of regression parameters. It is possible to pass a function
#' such as for example \code{sandwich::vcovHC} which will be called as:
#' \code{foo(object, ...)} and a user may specify additional arguments of a 
#' function in \code{...} argument. If not provided an estimate for covariance
#' matrix will be set by calling appropriate \code{vcov} method.
#' @param alpha Significance level for confidence intervals --
#' Either a single numeric value or a vector of length equal to number of 
#' sub populations specified in \code{stratas}. 
#' If missing it is set to \code{.05} in \code{singleR} method.
#' @param ... A vector of arguments to be passed to other functions.
#' For \code{singleR} method for this functions arguments in \code{...} are 
#' passed to either \code{cov} if argument provided was a function or 
#' \code{vcov} if \code{cov} argument was missing at call.
#' 
#' \loadmathjax
#' @details In single source capture-recapture models the most frequently used
#' estimate for population size is Horwitz-Thompson type estimate
#' 
#' \mjdeqn{\hat{N} = \sum_{k=1}^{N}\frac{I_{k}}{\mathbb{P}(Y_{k}>0)} = \sum_{k=1}^{N_{obs}}\frac{1}{1-\mathbb{P}(Y_{k}=0)}}{N = Sum_k=1^N I_k/P(Y_k > 0) = Sum_k=1^N_obs 1/(1-P(Y_k = 0))}
#'
#' where \mjeqn{I_{k}=I_{Y_{k} > 0}}{I_k=I_(Y_k > 0)} are indicator variables, 
#' with value 1 if kth unit was observed at least once and 0 otherwise and
#' the inverse probabilistic weights weights for units observed in the data
#' \mjeqn{\tfrac{1}{\mathbb{P}(Y_{k}>0)}}{1/P(Y_k > 0)}
#' are estimated using fitted linear predictors.
#' 
#' The estimates for different sub populations are made by changing the
#' \mjeqn{I_{k}=I_{Y_{k} > 0}}{I_k=I_(Y_k > 0)} indicator variables
#' to refer not to the population as a whole but to the sub populations that
#' are being considered i.e. by changing values from 1 to 0 if kth unit is not
#' a member of sub population that is being considered at the moment.
#' 
#' The estimation of variance for these estimates and estimation of variance for
#' estimate of population size for the whole population follow the same relation
#' as the one described above.
#' 
#' @seealso [vcov.singleR()] [estimatePopsize()]
#'
#' @return A \code{data.frame} object with row names being the names of specified 
#' sub populations either provided or inferred.
#' @export
stratifyPopEst <- function(object, stratas, alpha, ...) {
  UseMethod("stratifyPopEst")
}
#' @title Statistical tests of goodness of fit.
#'
#' @description Performs two statistical test on observed and fitted
#' marginal frequencies. For G test the test statistic is computed as:
#' \loadmathjax
#' \mjdeqn{G = 2\sum_{k}O_{k}\ln{\left(\frac{O_{k}}{E_{k}}\right)}}{G = 2 * Sum O_k ln (O_k/E_k)}
#' and for \mjeqn{\chi^{2}}{X2} the test statistic is computed as:
#' \mjdeqn{\chi^{2} = \sum_{k}\frac{\left(O_{k}-E_{k}\right)^{2}}{E_{k}}}{X2 = Sum (O_k - E_k)^2/E_k}
#' where \mjeqn{O_{k},E_{k}}{O_k,E_k} denoted observed and fitted frequencies respectively.
#' Both of these statistics converge to \mjeqn{\chi^2}{X2} distribution asymptotically 
#' with the same degrees of freedom.
#' 
#' The convergence of \mjeqn{G, \chi^2}{G, X2} statistics to \mjeqn{\chi^2}{X2}
#' distribution may be violated if expected counts in cells are too low, 
#' say < 5, so it is customary to either censor or omit these cells.
#' 
#' @param object Object of singleRmargin class.
#' @param df Degrees of freedom if not provided the function will try and manually
#' but it is not always possible.
#' @param dropl5 A character indicating treatment of cells with frequencies < 5 
#' either grouping them, droping or leaving them as is. Defaults to drop.
#' @param ... Currently does nothing.
#'
#' @method summary singleRmargin
#' @return A chi squared test and G test for comparison between fitted and observed marginal frequencies.
#' @examples 
#' # Create a simple model
#' Model <- estimatePopsize(formula = capture ~ ., 
#' data = netherlandsimmigrant, 
#' model = ztpoisson, 
#' method = "IRLS")
#' plot(Model, "rootogram")
#' # We see a considerable lack of fit
#' summary(marginalFreq(Model), df = 1, dropl5 = "group")
#' @exportS3Method
summary.singleRmargin <- function(object, df,
                                  dropl5 = c("drop", 
                                             "group", 
                                             "no"), 
                                  ...) {
  if (missing(dropl5)) {dropl5 <- "drop"}
  y <- object$y
  if (grepl("zot", object$name) & (1 %in% names(y))) {y <- y[-1]}
  A <- object$table[names(y)]
  if ((missing(df)) && (object$df < 1)) {
    warning("Degrees of freedom may be inacurate.")
    df <- 1
  } else if (missing(df)) {
    df <- object$df
  }
  if(dropl5 == "group") {
    l <- (A < 5)
    if(!all(l == FALSE)) {
      y <- c(y[!l], sum(y[l]))
      A <- c(A[!l], sum(A[l]))
    }
  } else if(dropl5 == "drop") {
    l <- (A < 5)
    y <- y[!l]
    A <- A[!l]
  }
  
  X2 <- sum(((A - y) ^ 2) / A)
  G <- 2 * sum(y * log(y / A))
  pval <- stats::pchisq(q = c(X2, G), df = df, lower.tail = FALSE)
  vect <- data.frame(round(c(X2, G), digits = 2),
                     rep(df, 2), signif(pval, digits = 2))
  rownames(vect) <- c("Chi-squared test", "G-test")
  colnames(vect) <- c("Test statistics", "df", "P(>X^2)")
  structure(
    list(Test = vect,
         l5 = switch(dropl5,
                     drop = "dropped",
                     group = "grouped",
                     no = "preserved"),
         y = y),
    class = "summarysingleRmargin"
  )
}
#' @title Obtain Covariance Matrix estimation.
#' 
#' @description A \code{vcov} method for \code{singleR} class.
#' 
#' @param object Object of singleR class.
#' @param type Type of estimate for covariance matrix for now either
#' expected (Fisher) information matrix or observed information matrix.
#' @param ... Additional optional arguments passed to following functions:
#' \itemize{
#' \item \code{solve} -- for inverting information matrixes.
#' \item \code{model.frame.singleR} -- forcreation of Xvlm matrix.
#' }
#' @details  Returns a estimated covariance matrix for model coefficients
#' calculated from analytic hessian or Fisher information matrix. 
#' Covariance type is taken from control parameter that have been provided
#' on call that created \code{object} if argumetns \code{type} was not specified.
#' 
#' @method vcov singleR
#' @return A covariance matrix for fitted coefficients, rows and columns of which 
#' correspond to parameters returned by \code{coef} method.
#' @seealso [vcovHC.singleR()] [sandwich::sandwich()]
#' @exportS3Method
vcov.singleR <- function(object, 
                         type = c("Fisher", 
                                  "observedInform"), 
                         ...) {
  if (missing(type) ){type <- object$populationSize$control$covType}
  X <- model.frame.singleR(object, ...)
  X <- subset(X, select = colnames(X)[-(attr(object$terms, "response"))], subset = object$which$reg)
  X <- singleRinternalGetXvlmMatrix(
    X = X, 
    nPar = object$model$parNum, 
    formulas = object$formula, 
    parNames = object$model$etaNames
  )
  res <- switch(
    type,
    "observedInform" = solve(
      -object$model$makeMinusLogLike(y = object$y[object$which$reg], X = X,
      weight = object$priorWeights[object$which$reg], deriv = 2)(object$coefficients),
      ...
    ),
    "Fisher" = {
      if (object$call$method == "IRLS") {W <- object$weights} else {W <- object$model$Wfun(prior = object$priorWeights[object$which$reg], eta = object$linearPredictors)};
      solve(
      singleRinternalMultiplyWeight(X = X, W = W) %*% X,
      ...
    )}
  )
  dimnames(res) <- list(names(object$coefficients), names(object$coefficients))
  res
}
#' @title Diagonal elements of the projection matrix
#' 
#' @description A method for \code{singleR} class for extracting diagonal 
#' elementrs of projection matrix.
#' 
#' @param model Object of singleR class.
#' @param ... Currently does nothing.
#' 
#' \loadmathjax
#' @details Since \code{singleRcapture} contains not only regular glm's but also
#' vglm's the hatvalues returns a matrix with number of columns corresponding to
#' number of linear predictors in a model, where kth column corresponds
#' to elements of the diagonal of projection matrix associated with kth linear 
#' predictor. For glm's  
#' \mjdeqn{\boldsymbol{W}^{\frac{1}{2}}\boldsymbol{X}\left(\boldsymbol{X}^{T}\boldsymbol{W}\boldsymbol{X}\right)^{-1}\boldsymbol{X}^{T}\boldsymbol{W}^{\frac{1}{2}}}{sqrt(W)X(X'WX)^-1X'sqrt(W)}
#' where \mjeqn{\boldsymbol{W}=\mathbb{E}\left(\text{Diag}\left(\frac{\partial^{2}\ell}{\partial\boldsymbol{\eta}^{T}\partial\boldsymbol{\eta}}\right)\right)}{W = E(diag(d^2.ln(L)/d.eta^2))} and \mjeqn{\boldsymbol{X}}{X} is a model (lm) matrix. 
#' For vglm's it is instead :
#' \mjdeqn{\boldsymbol{X}_{vlm}\left(\boldsymbol{X}_{vlm}^{T}\boldsymbol{W}\boldsymbol{X}_{vlm}\right)^{-1}\boldsymbol{X}_{vlm}^{T}\boldsymbol{W}}{X_vlm(X_vlm'WX_vlm)^-1X_vlm'W}
#' where 
#' \mjdeqn{
#' \boldsymbol{W} = \mathbb{E}\left(\begin{bmatrix}
#' \text{Diag}(\frac{\partial^{2}\ell}{\partial\eta_{1}^{T}\partial\eta_{1}}) &
#' \text{Diag}(\frac{\partial^{2}\ell}{\partial\eta_{1}^{T}\partial\eta_{2}}) &
#' \dotso & \text{Diag}(\frac{\partial^{2}\ell}{\partial\eta_{1}^{T}\partial\eta_{p}})\cr
#' \text{Diag}(\frac{\partial^{2}\ell}{\partial\eta_{2}^{T}\partial\eta_{1}}) &
#' \text{Diag}(\frac{\partial^{2}\ell}{\partial\eta_{2}^{T}\partial\eta_{2}}) &
#' \dotso & \text{Diag}(\frac{\partial^{2}\ell}{\partial\eta_{2}^{T}\partial\eta_{p}})\cr
#' \vdots & \vdots & \ddots & \vdots\cr
#' \text{Diag}(\frac{\partial^{2}\ell}{\partial\eta_{p}^{T}\partial\eta_{1}}) &
#' \text{Diag}(\frac{\partial^{2}\ell}{\partial\eta_{p}^{T}\partial\eta_{2}}) &
#' \dotso & \text{Diag}(\frac{\partial^{2}\ell}{\partial\eta_{p}^{T}\partial\eta_{p}})
#' \end{bmatrix}\right)}{W = E(matrix(
#' diag(d^2.ln(L)/d.eta_1^2), diag(d^2.ln(L)/d.eta_1 d.eta_2), ..., diag(d^2.ln(L)/d.eta_1 d.eta_p)
#' diag(d^2.ln(L)/d.eta_2 d.eta_1), diag(d^2.ln(L)/d.eta_2^2), ..., diag(d^2.ln(L)/d.eta_2 d.eta_p)
#' .............................................................................................
#' diag(d^2.ln(L)/d.eta_p d.eta_1), diag(d^2.ln(L)/d.eta_p d.eta_2), ..., diag(d^2.ln(L)/d.eta_p^2)))}
#' is a block matrix constructed by taking the expected  value from diagonal 
#' matrixes corresponding to second derivatives with respect to each linear 
#' predictor (and mixed derivatives) and 
#' \mjeqn{\boldsymbol{X}_{vlm}}{X_vlm} is a model (vlm) matrix constructed using 
#' "main" formula and additional formulas specified in \code{controlModel}. 
#' @method hatvalues singleR
#' @importFrom stats hatvalues
#' @return A matrix with n rows and p columns where n is a number of observations
#' in the data and k is number of regression parameters.
#' @exportS3Method 
hatvalues.singleR <- function(model, ...) {
  X <- model.frame.singleR(model, ...)
  X <- subset(X, select = colnames(X)[-(attr(model$terms, "response"))], subset = model$which$reg)
  X <- singleRinternalGetXvlmMatrix(X = X, nPar = model$model$parNum, formulas = model$formula, parNames = model$model$etaNames)
  if (isTRUE(model$call$method == "IRLS")) {
    W <- model$weights
  } else {
    W <- model$model$Wfun(prior = model$priorWeights[model$which$reg], eta = if (model$model$family == "zelterman") model$linearPredictors[model$which$reg, ] else model$linearPredictors)
  }
  mlt <- singleRinternalMultiplyWeight(X = X, W = W)
  hatvalues <- diag(X %*% solve(mlt %*% X) %*% mlt)
  hatvalues <- matrix(hatvalues, ncol = model$model$parNum, dimnames = list(1:(length(hatvalues) / model$model$parNum), model$model$etaNames))
  hatvalues
}
#' @title Regression deletion diagnostics.
#' \loadmathjax
#' @description A \code{dfbeta} method for \code{singleR} class.
#' 
#' @param model Fitted object of singleR class
#' @param maxitNew Maximal number of iterations for regressions with starting points
#' \mjeqn{\hat{\boldsymbol{\beta}}}{beta} on data specified at call for \code{model}
#' after the romoval of kth row. By default 1.
#' @param ... Additional optional arguments passed to the following functions:
#' \itemize{
#' \item \code{controlMethod} -- For controlling the simulation of dfbeta.
#' }
#' 
#' @return A matrix with n rows and p observations where p is a number of
#' units in data and p is the number of regression parameters.
#' 
#' K'th row of this matrix corresponds to 
#' \mjeqn{\hat{\boldsymbol{\beta}}-\hat{\boldsymbol{\beta}}_{-k}}{beta-beta_(-k)}
#' where \mjeqn{\hat{\boldsymbol{\beta}}_{-k}}{beta_(-k)} is a vector of estimates
#' for regression parameters after the removal of k'th row from the data.
dfbetasingleR <- function(model,
                          maxitNew = 1,
                          ...) {
  # formula method removed since it doesn't give good results will reimplement if we find better formula
  X <- model.frame.singleR(model, ...)
  y <- if (is.null(model$y)) stats::model.response(X) else model$y
  X <- subset(X, select = colnames(X)[-(attr(model$terms, "response"))], subset = model$which$reg)
  y <- y[model$which$reg]
  cf <- model$coefficients
  pw <- model$priorWeights[model$which$reg]
  res <- matrix(nrow = nrow(X), ncol = length(cf))
  for (k in 1:nrow(X)) {
    res[k, ] <- cf - estimatePopsize.fit(
      control = controlMethod(
        silent = TRUE, 
        start = cf,
        maxiter = maxitNew + 1,
        ...
      ),
      y = y[-k],
      X = singleRinternalGetXvlmMatrix(X = subset(X, rownames(X) != rownames(X)[k]), nPar = model$model$parNum, formulas = model$formula, parNames = model$model$etaNames),
      start = cf,
      family = model$model,
      priorWeights = pw[-k],
      method = model$call$method
    )$beta
  }
  colnames(res) <- names(model$coefficients)
  res
}
#' @title Confidence Intervals for Model Parameters
#' 
#' @description A function that computes studentized confidence intervals
#' for model coefficients.
#' 
#' @param object Object of singleR class.
#' @param parm Names of parameters for which confidence intervals are to be 
#' computed, if missing all parameters will be considered.
#' @param level Confidence level for intervals.
#' @param ... Currently does nothing.
#' 
#' @method confint singleR
#' @return An object with named columns that include upper and 
#' lower limit of confidence intervals.
#' @exportS3Method
confint.singleR <- function(object,
                            parm, 
                            level = 0.95, 
                            ...) {
  std <- sqrt(diag(vcov.singleR(object)))
  if (missing(parm)) {
    coef <- object$coefficients
  } else {
    coef <- object$coefficients[parm]
    std <- std[parm]
  }
  sc <- qnorm(p = 1 - (1 - level) / 2)
  res <- data.frame(coef - sc * std, coef + sc * std)
  colnames(res) <- c(paste0(100 * (1 - level) / 2, "%"),
                     paste0(100 * (1 - (1 - level) / 2), "%"))
  res
}

# There is no need for doccumenting the following methods:
#' @method residuals singleR
#' @importFrom stats residuals
#' @exportS3Method
residuals.singleR <- function(object,
                              type = c("pearson",
                                       "pearsonSTD",
                                       "response",
                                       "working",
                                       "deviance",
                                       "all"),
                              ...) {
  type <- match.arg(type)
  res <- object$residuals
  wts <- object$priorWeights
  mu <- object$fitt.values
  y <- if (is.null(object$y)) stats::model.response(model.frame(object)) else object$y
  if (!(all(object$which$reg == object$which$est)) && type == "all") stop("type = all is not aviable for some models")
  if (type == "pearsonSTD" && object$model$parNum > 1) {stop("Standardized pearson residuals not yet implemented for models with multiple linear predictors")}
  rs <- switch(
    type,
    working = as.data.frame(object$model$funcZ(eta = if (object$model$family == "zelterman") object$linearPredictors[object$which$reg, ] else object$linearPredictors, weight = object$weights, y = y[object$which$reg]), col.names = paste0("working:", object$model$etaNames)),
    response = res,
    pearson = data.frame("pearson" = (if (object$model$family == "zelterman") res$mu[object$which$reg] else res$mu) / sqrt(object$model$variance(eta = if (object$model$family == "zelterman") object$linearPredictors[object$which$reg, ] else object$linearPredictors, type = "trunc"))),
    pearsonSTD = data.frame("pearsonSTD" = (if (object$model$family == "zelterman") res$mu[object$which$reg] else res$mu) / sqrt((1 - hatvalues(object)) * object$model$variance(eta = if (object$model$family == "zelterman") object$linearPredictors[object$which$reg, ] else object$linearPredictors, type = "trunc"))),
    deviance = data.frame("deviance" = object$model$devResids(y = y[object$which$reg], eta = object$linearPredictors, wt = wts[object$which$reg])),
    all = {colnames(res) <- c("muResponse", "linkResponse");
      data.frame(
      as.data.frame(object$model$funcZ(eta = if (object$model$family == "zelterman") object$linearPredictors[object$which$reg, ] else object$linearPredictors, weight = object$weights, y = y[object$which$reg]), col.names = paste0("working:", object$model$etaNames)),
      res,
      "pearson" = as.numeric((if (object$model$family == "zelterman") res$mu[object$which$reg] else res$mu) / sqrt(object$model$variance(eta = if (object$model$family == "zelterman") object$linearPredictors[object$which$reg, ] else object$linearPredictors, type = "trunc"))),
      "pearsonSTD" = if (object$model$parNum == 1) as.numeric((if (object$model$family == "zelterman") res$mu[object$which$reg] else res$mu) / sqrt((1 - hatvalues(object)) * object$model$variance(eta = if (object$model$family == "zelterman") object$linearPredictors[object$which$reg, ] else object$linearPredictors, type = "trunc"))) else 0,
      "deviance" = as.numeric(object$model$devResids(y = y[object$which$reg], eta = object$linearPredictors, wt = wts[object$which$reg])),
      row.names = rownames(object$linearPredictors)
    )}
  )
  rs
}
#' @method family singleR
#' @importFrom stats family
#' @exportS3Method
family.singleR <- function(object, ...) {
  object$model
}
#' @method print summarysingleRmargin
#' @exportS3Method 
print.summarysingleRmargin <- function(x, ...) {
  cat("Test for Goodness of fit of a regression model:\n",
      "\n", sep = "")
  print(x$Test)
  cat("\n--------------------------------------------------------------",
      "\nCells with fitted frequencies of < 5 have been", x$l5, 
      "\nNames of cells used in calculating test(s) statistic:", names(x$y), 
      "\n", sep = " ")
}
#' @method AIC singleR
#' @importFrom stats AIC
#' @exportS3Method 
AIC.singleR <- function(object, ...) {
  2 * (length(object$coefficients) - object$logL)
}
#' @method BIC singleR
#' @importFrom stats BIC
#' @exportS3Method 
BIC.singleR <- function(object, ...) {
  length(object$coefficients) * log(sum(object$which$reg)) - 2 * object$logL
}
#' @method extractAIC singleR
#' @importFrom stats extractAIC
#' @exportS3Method 
extractAIC.singleR <- function(fit, scale, k = 2, ...) {
  -2 * fit$logL + k * length(fit$coefficients)
}
#' @method dfbeta singleR
#' @importFrom stats dfbeta
#' @exportS3Method 
dfbeta.singleR <- function(model, ...) {
  dfbetasingleR(model, ...)
}
#' @method logLik singleR
#' @importFrom stats logLik
#' @exportS3Method 
logLik.singleR <- function(object, ...) {
  val <- object$logL
  attr(val, "nobs") <- dim(residuals(object))[1]
  attr(val, "df") <- length(object$coefficients)
  class(val) <- "logLik"
  val
}
#' @method model.frame singleR
#' @importFrom stats glm
#' @importFrom stats model.frame
#' @exportS3Method 
model.frame.singleR <- function(formula, ...) {
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), 
                      names(dots), 
                      0L)]
  if (length(nargs) || is.null(formula$modelFrame)) {
    fcall <- formula$call
    fcall$method <- "model.frame"
    fcall[[1L]] <- quote(stats::glm)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    eval(fcall, env)
  }
  else formula$modelFrame
}

#' @method model.matrix singleR
#' @importFrom stats model.matrix
#' @exportS3Method 
model.matrix.singleR <- function(object, type = c("lm", "vlm"), ...) {
  if (missing(type)) type <- "lm"
  switch (type,
    lm = {
      if (is.null(object$X)) {
        X <- model.frame(object);
        X <- model.matrix(object$terms, X)
      } else {
        X <- object$X
      }
      subset(X, subset = object$which$reg)
      },
    vlm = {
      X <- model.frame.singleR(object, ...);
      X <- subset(X, select = colnames(X)[-(attr(object$terms, "response"))], subset = object$which$reg);
      singleRinternalGetXvlmMatrix(X = X, nPar = object$model$parNum, 
      formulas = object$formula, parNames = object$model$etaNames);
      }
  )
}
#' @method redoPopEstimation singleR
#' @rdname redoPopEstimation
#' @exportS3Method
redoPopEstimation.singleR <- function(object, cov = NULL, ...) {
  Xvlm <- model.matrix(object, "vlm")
  singleRcaptureinternalpopulationEstimate(
    y = object$y[object$which$reg],
    formulas = object$formula,
    X = model.matrix(object),
    grad = object$model$makeMinusLogLike(y = object$y[object$which$reg], 
    X = Xvlm, weight = object$priorWeights, deriv = 1),
    hessian = object$model$makeMinusLogLike(y = object$y[object$which$reg], 
    X = Xvlm, weight = object$priorWeights, deriv = 2),
    popVar = if (is.null(object$call$popVar)) "analytic" else object$call$popVar,
    weights = object$priorWeights[object$which$reg],
    eta = object$linearPredictors,
    family = object$model,
    beta = object$coefficients,
    control = object$populationSize$control,
    Xvlm = Xvlm,
    W = if (object$call$method == "IRLS") object$weights else object$model$Wfun(prior = object$priorWeights, eta = object$linearPredictors),
    sizeObserved = object$sizeObserved,
    modelFrame = model.frame.singleR(object, ...),
    cov = cov
  )
}
#' @method dfpopsize singleR
#' @rdname dfpopsize
#' @exportS3Method 
dfpopsize.singleR <- function(model, dfbeta = NULL, observedPop = FALSE, ...) {
  if (isTRUE(model$call$popVar == "bootstrap")) warning("dfpopsize may (in some cases) not work correctly when bootstrap was chosen as population variance estimate.")
  dfb <- if (is.null(dfbeta)) {dfbeta(model, ...)} else {dfbeta}
  if (model$model$family == "zelterman") {
    dfbnew <- matrix(0, ncol = ncol(dfb), nrow = NROW(model$priorWeights))
    rownames(dfbnew) <- as.character(1:NROW(dfbnew))
    dfbnew[model$which$reg, ] <- dfb
    dfb <- dfbnew
  }
  X <- model.frame.singleR(model, ...)
  X <- subset(X, select = colnames(X)[-(attr(model$terms, "response"))], subset = model$which$est)
  N <- model$populationSize$pointEstimate
  range <- 1:NROW(dfb)
  res <- vector("numeric", length = NROW(dfb))
  pw <- model$priorWeights[model$which$est]
  for (k in range) {
    cf <- model$coefficients - dfb[k, ]
    res[k] <- model$model$pointEst(
      eta = matrix(singleRinternalGetXvlmMatrix(X = subset(X, rownames(X) != rownames(X)[k]), nPar = model$model$parNum, formulas = model$formula, parNames = model$model$etaNames) %*% cf, ncol = model$model$parNum),
      pw = pw[-k]) + model$trcount
  }
  
  if(isTRUE(observedPop) & (grepl("zot", model$model$family) | model$model$family == "chao")) {
    res1 <- vector(mode = "numeric", length = model$sizeObserved)
    names(res1) <- 1:model$sizeObserved
    res1[model$which$est] <- N - res
    res1[!model$which$est] <- 1
    
  } else {
    res1 <- N - res
  }
  
  res1
}
#' @method stratifyPopEst singleR
#' @rdname stratifyPopEst
#' @importFrom stats vcov
#' @exportS3Method
stratifyPopEst.singleR <- function(object, stratas, alpha, cov = NULL, ...) {
  # if stratas is unspecified get all levels of factors in modelFrame
  if (missing(stratas)) {
    stratas <- names(which(attr(object$terms, "dataClasses") == "factor"))
    stratas <- stratas[stratas %in% attr(object$terms, "term.labels")]
  }
  # if significance level is unspecified set it to 5%
  if (missing(alpha)) alpha <- .05
  
  # convert stratas to list for all viable types of specifying the argument
  if (inherits(stratas, "formula")) {
    modelFrame <- model.frame(object)
    # TODO
  } else if (is.list(stratas)) {
    if (!all(sapply(stratas, is.logical))) {
      stop("Invalid way of specifying subpopulations in stratas. If stratas argument is a list ")
    }
    if (length(stratas[[1]]) != object$sizeObserved) stop("Elements of stratas object should have length equal to number of observed units.")
  } else if (is.logical(stratas)) {
    if (length(stratas) != object$sizeObserved) stop("Stratas object should have length equal to number of observed units.")
    stratas <- list(strata = stratas)
  } else if (is.character(stratas)) {
    modelFrame <- model.frame(object)
    out <- list()
    for (k in stratas) {
      if (!(k %in% colnames(modelFrame))) stop("Variable specified in stratas is not present in model frame.")
      if (!(is.factor(modelFrame[,k])) & !(is.character(modelFrame[,k]))) stop("Variable specified in stratas is not a factor or a character vector.")

      if (is.factor(modelFrame[, k])) {
        # this makes a difference on factor that is not present
        for (t in levels(modelFrame[, k])) {
          out[[paste0(as.character(k), "==", t)]] <- (modelFrame[, k] == t)
        }
      } else {
        for (t in unique(modelFrame[, k])) {
          out[[paste0(as.character(k), "==", t)]] <- (modelFrame[, k] == t)
        }
      }
    }
    stratas <- out
  } else {
    # a formula, a list with logical vectors specifying different sub populations or a single logical vector or a vector with names of factor variables.
    errorMessage <- paste0("Invalid way of specifying subpopulations in stratas.\n", 
    "Please provide either:\n",
    "(1) - a list with logical vectors specifying different sub populations\n",
    "(2) - a single logical vector\n",
    "(3) - a formula\n",
    "(4) - a vector with names of factor variables\n",)
    stop(errorMessage)
  }
  
  # get neccesary model info AFTER possible error in function
  family <- family(object = object, ...)
  priorWeights <- object$priorWeights
  eta <- object$linearPredictors
  Xvlm <- model.matrix(object, "vlm")
  
  # get covariance matrix
  if (is.function(cov)) cov <- cov(object, ...)
  if (is.null(cov)) cov <- vcov(object, ...)
  
  obs <- vector(mode = "numeric", length = length(stratas))
  est <- vector(mode = "numeric", length = length(stratas))
  stdErr <- vector(mode = "numeric", length = length(stratas))
  cnfStudent <- matrix(nrow = length(stratas), ncol = 2)
  cnfChao <- matrix(nrow = length(stratas), ncol = 2)
  sc <- qnorm(p = 1 - alpha / 2)
  if (length(sc) != length(stratas)) sc <- rep(sc, length.out = length())
  
  #TODO adjust for chao
  for (k in 1:length(stratas)) {
    cond <- stratas[[k]]
    obs[k] <- sum(cond)
    if (obs[k] > 0) {
      est[k] <- family$pointEst(pw = priorWeights[cond], eta = eta[cond, ])
      stdErr[k] <- family$popVar(pw = priorWeights[cond], eta = eta[cond, ], cov = cov, Xvlm = Xvlm[rep(cond, length(family$etaNames)), ]) ^ .5
      cnfStudent[k, ] <- est[k] + c(-sc * stdErr[k], sc * stdErr[k])
      G <- exp(sc * sqrt(log(1 + (stdErr[k]^2) / ((est[k] - obs[k]) ^ 2))))
      cnfChao[k, ] <- obs[k] + c((est[k] - obs[k]) / G, (est[k] - obs[k]) * G)
    } else {
      est[k] <- 0
      stdErr[k] <- 0
      cnfStudent[k, ] <- c(0, 0)
      cnfChao[k, ] <- c(0, 0)
    }
  }
  
  result <- data.frame(
    obs, est, 100 * obs / est, stdErr, 
    cnfStudent[, 1], cnfStudent[, 2], 
    cnfChao[, 1], cnfChao[, 2],
    row.names = names(stratas)
  )
  if (length(unique(alpha)) == 1) {
    nma <- c(as.character(100 * alpha / 2), as.character(100 * (1 - alpha / 2)))
  } else {
    nma <- c("LowerBounds", "UpperBounds")
  }
  colnames(result) <- c("Observed", "Estimated", "ObservedProcentage", "StdError",
  paste0("Studentized - ", nma, "%"), paste0("Chao - ", nma, "%"))
  
  result
}
#' @method popSizeEst singleR
#' @rdname popSizeEst
#' @exportS3Method
popSizeEst.singleR <- function(object, ...) {
  object$populationSize
}
#' @importFrom stats cooks.distance
#' @method cooks.distance singleR
#' @exportS3Method 
cooks.distance.singleR <- function(model, ...) {
  if (model$model$parNum > 1) stop("Cooks distance is only implemented for single parameter families.")
  res <- residuals(model, type = "pearsonSTD") ^ 2
  res <- res[, 1]
  ht <- hatvalues(model)
  res <- (res * (ht / (length(model$coefficients))))
  names(res) <- rownames(ht)
  res
}
#' @method print popSizeEstResults
#' @exportS3Method 
print.popSizeEstResults <- function(x, ...) {
  cat("Point estimate: ", x$pointEstimate, "\nVariance: ", x$variance, "\n", (1-x$control$alpha) * 100, "% confidence intervals:\n", sep = "")
  print(x$confidenceInterval)
  invisible(x)
}
#' @method print summarysingleR
#' @importFrom stats printCoefmat
#' @exportS3Method 
print.summarysingleR <- function(x, 
                                 signif.stars = getOption("show.signif.stars"), 
                                 digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Pearson Residuals:\n")
  print(summary(c(x$residuals[, 1])))
  cat("\nCoefficients:\n")
  #print(ob)
  cond <- sapply(x$model$etaNames, 
  FUN = function(k) {sapply(strsplit(rownames(x$coefficients), split = ":"), 
                     FUN = function(x) {x[[length(x)]] == k})})
  for (k in x$model$etaNames) {
    if (!all(cond[,k] == FALSE)) {
      cat("-----------------------\nFor linear predictors associated with:", k, "\n")
      # Base liblary has no str_length and no proper sub_string this is a slightly convoluted way of making do without them
      toPrint <- subset(x$coefficients, cond[,k])
      lengths <- sapply(rownames(toPrint), function(x) {length(strsplit(x, split = "")[[1]])})
      lK <- length(unlist(strsplit(k, split = "")))
      rownames(toPrint) <- sapply(1:nrow(toPrint), function(x) {substr(x = rownames(toPrint)[x], start = 1, stop = lengths[[x]] - (1 + lK))})
      #print(toPrint)
      printCoefmat(toPrint, digits = digits, 
                   signif.stars = signif.stars, 
                   signif.legend = if (k == x$model$etaNames[length(x$model$etaNames)]) signif.stars else FALSE, 
                   P.values = TRUE, has.Pvalue = TRUE, 
                   na.print = "NA", ...)
    } else {
      cat("-----------------------\nFor linear predictors associated with:", k, "\n")
      #print(subset(x$coefficients, rowSums(cond) == 0))
      printCoefmat(subset(x$coefficients, rowSums(cond) == 0), 
                   digits = digits, 
                   signif.stars = signif.stars, 
                   signif.legend = if (k == x$model$etaNames[length(x$model$etaNames)]) signif.stars else FALSE, 
                   P.values = TRUE, has.Pvalue = TRUE, 
                   na.print = "NA", ...)
    }
  }
  
  cat("\n")
  if (!is.null(x$correlation)) {
    corr <- round(x$correlation, digits = 2)
    corr[!lower.tri(corr)] <- ""
    print(corr[-1, -dim(corr)[2]], quote = FALSE)
    cat("\n")
  }
  sd <- sqrt(x$populationSize$variance)
  if (x$populationSize$control$sd == "normalMVUE") {
    sd <- sd / (sqrt(2 / (x$sizeObserved - 1)) * exp(lgamma(x$sizeObserved / 2) - lgamma((x$sizeObserved - 1) / 2)))
  }
  cat("AIC: ", x$aic,
      "\nBIC: ", x$bic,
      "\nResidual deviance: ", x$deviance,
      "\n\nLog-likelihood: ", x$logL, " on ", x$dfResidual, " Degrees of freedom ",
      if (isTRUE(x$call$method == "IRLS")) {
        "\nNumber of iterations: "
      } else {
        "\nNumber of calls to log-likelihood function: " # optim does not allow for accessing information
        # on number of iterations performed only a number of calls for gradient and objective function
      }, x$iter[1], 
      "\n-----------------------",
      "\nPopulation size estimation results: ",
      "\nPoint estimate ", x$populationSize$pointEstimate, 
      "\nObserved proportion: ", round(100 * x$sizeObserved / x$populationSize$pointEstimate, digits = 1), "% (N obs = ", x$sizeObserved, ")",
      if (!is.null(x$skew)) {"\nBoostrap sample skewness: "}, if (!is.null(x$skew)) {x$skew}, if (!is.null(x$skew)) {"\n0 skewness is expected for normally distributed vairable\n---"},
      if (isTRUE(x$call$popVar == "bootstrap")) {"\nBootstrap Std. Error "} else {"\nStd. Error "}, sd,
      "\n", (1 - x$populationSize$control$alpha) * 100, "% CI for the population size:\n", sep = "")
  print(x$populationSize$confidenceInterval)
  cat((1 - x$populationSize$control$alpha) * 100, "% CI for the share of observed population:\n", sep = "")
  dd <- as.data.frame(x$populationSize$confidenceInterval)
  if (ncol(dd) == 1) {
    vctpop <- sort(x$populationSize$confidenceInterval, decreasing = TRUE)
    names(vctpop) <- rev(names(vctpop))
    print(100 * x$sizeObserved / vctpop)
  } else {
    print(data.frame(
      lowerBound = 100 * x$sizeObserved / dd[, 2], 
      upperBound = 100 * x$sizeObserved / dd[, 1],
      row.names = rownames(dd)
    ))
  }
  invisible(x)
}
#' @method print singleR
#' @exportS3Method 
print.singleR <- function(x, ...) {
  cat("Call: ")
  print(x$call)
  cat("\nCoefficients:\n")
  print(coef(x))
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e NULL);", x$dfResidual, "Residual")
  cat("\nAIC: ", signif(AIC(x)), "\nBIC: ", signif(BIC(x)), "\nResidual deviance: ", signif(x$deviance))
  cat("\n-----------------------------------\nPopulation size estimation results:\n")
  print(x$populationSize)
  invisible(x)
}

#' @importFrom stats fitted
#' @method fitted singleR
#' @exportS3Method 
fitted.singleR <- function(object,
                           type = c("mu", 
                                    "link"), # maybe add marginal frequencies/other options
                           ...) {
  if (missing(type)) type <- "mu"
  object$fitt.values[[type]] # fitted should return either E(Y) or E(Y|Y>0) otherwise we're breaking R conventions
}
# TODO:: update
#' #' simulate
#' #' 
#' #' An S3class for \code{stats::simulate} to handle \code{singleR} objects.
#' #' @author Maciej Beręsewicz
#' #' 
#' #' @param object an object representing a fitted model.
#' #' @param nsim number of response vectors to simulate. Defaults to \code{1}.
#' #' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).
#' #' @param ... additional optional arguments.
#' #' @return a \code{data.frame} with \code{n} rows and \code{nsim} columns.
#' #' @seealso [stats::simulate()]
#' #' @examples 
#' #' set.seed(1)
#' #' N <- 10000
#' #' gender <- rbinom(N, 1, 0.2)
#' #' eta <- -1 + 0.5*gender
#' #' counts <- rpois(N, lambda = exp(eta))
#' #' df <- data.frame(gender, eta, counts)
#' #' df2 <- subset(df, counts > 0)
#' #' mod1 <-  estimatePopsize(formula = counts ~ 1 + gender, data = df2, 
#' #' model = "ztpoisson", method = "optim", popVar = "analytic")
#' #' mod1_sims <- simulate(mod1, nsim=10)
#' #' colMeans(mod1_sims) 
#' #' @importFrom stats simulate
#' #' @method simulate singleR
#' #' @exportS3Method
#' simulate.singleR <- function(object, nsim=1, seed = NULL, ...) {
#'   ###################
#'   ### TODO update ###
#'   ###################
#'   stop("Not yet updated")
#'   if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
#'     runif(1)
#'   if (is.null(seed))
#'     RNGstate <- get(".Random.seed", envir = .GlobalEnv)
#'   else {
#'     R.seed <- get(".Random.seed", envir = .GlobalEnv)
#'     set.seed(seed)
#'     RNGstate <- structure(seed, kind = as.list(RNGkind()))
#'     on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
#'   }
#'   
#'   linpred <- object$linearPredictors
#'   n <- NROW(linpred)
#'   if (grepl("negbin", object$model$family)) {
#'     val <- object$model$simulate(n*nsim, exp(linpred), exp(-object$dispersion)) # dispersion argument will be removed in next update remember to change that to accomodate that
#'   } else {
#'     val <- object$model$simulate(n*nsim, exp(linpred))
#'   }
#'   
#'   dim(val) <- c(n, nsim)
#'   val <- as.data.frame(val)
#'   names(val) <- paste0("sim_", seq_len(nsim))
#'   return(val)
#' }
