#' @title Observed and fitted marginal Frequencies
#' @author Piotr Chlebicki
#' \loadmathjax
#' @description A function that given a fitted \code{singleR} class object 
#' computed marginal frequencies by as sum of probability density functions
#' for each unit in data at each point i.e. kth element of marginal frequency
#' table is given by \mjseqn{\sum_{j=1}^{N_{obs}}\mathbb{P}(Y_{j}=k|\eta_{j})}. 
#' For k=0 only (if specified at call) they are computed as 
#' \mjseqn{\hat{N}-N_{obs}} because 
#' \mjseqn{\boldsymbol{f}_{0}} is assumed to 
#' the unobserved part of the studied population.
#' 
#' These frequencies are useful in diagnostics for count data regression, such
#' as assessment of fit.
#'
#' @param object object of \code{singleR} class.
#' @param includeones logical value indicating whether to include the estimated number of zero counts.
#' @param includezeros logical value indicating whether to include one counts in the zero-one truncated models.
#' @param onecount a numeric value indicating number of one counts if null \code{trcount} from object will be assumed to be a number one counts.
#' @param range optional argument specifying range of selected Y values.
#' @param ... currently does nothing.
#'
#' @return A list with observed name of the fitted model family degrees of freedom and observed and fitted marginal frequencies.
#' @seealso [estimatePopsize()] -- where example of usage is provided
#' @export
marginalFreq <- function(object,
                         includeones = TRUE, # matters only for zero one truncated models
                         includezeros = TRUE,
                         onecount = NULL,
                         range,
                         ...) {
  y <- if (is.null(object$y)) stats::model.response(model.frame(object)) else object$y
  if (missing(range)) {range <- (min(y):max(y))}
  y <- table(y)[names(table(y)) %in% as.character(range)]
  y <- y[!is.na(y)]
  
  trcount <- object$trcount
  if(!is.null(onecount)) {
    trcount <- onecount
  }
  
  # PMF for truncated distributions:
  probFun <- object$model$densityFunction
  
  res <- sapply(
    1:nrow(object$linearPredictors), 
    FUN = function(x) {object$model$densityFunction(
      x = range, 
      eta = object$linearPredictors[x, , drop = FALSE]
    )}
  )
  res <- rowSums(res)
  names(res) <- as.character(range)
  
  if(isTRUE(includeones) & grepl(object$model$family, pattern = "^zot")) {
    res <- c(trcount, res)
    names(res)[1] <- "1"
    y <- c(trcount, y)
    names(y)[1] <- "1"
    if (!is.null(onecount)) {res[1] <- onecount; y[1] <- onecount}
  }
  
  if(isTRUE(includezeros)) {
    res <- c(object$populationSize$pointEstimate - length(object$y), res)
    names(res)[1] <- "0"
  }
  res <- structure(list(
    table = res, y = y, 
    df = length(y) - length(object$coefficients) - 1,
    name = object$model$family
  ),
  class = c("singleRmargin"))
  res
}

# moved from methods. 
# TODO:: add simulated p-value
#' @title Statistical tests of goodness of fit.
#'
#' @description Performs two statistical test on observed and fitted
#' marginal frequencies. For G test the test statistic is computed as:
#' \loadmathjax
#' \mjsdeqn{G = 2\sum_{k}O_{k}\ln{\left(\frac{O_{k}}{E_{k}}\right)}}
#' and for \mjseqn{\chi^{2}} the test statistic is computed as:
#' \mjsdeqn{\chi^{2} = \sum_{k}\frac{\left(O_{k}-E_{k}\right)^{2}}{E_{k}}}
#' where \mjseqn{O_{k},E_{k}} denoted observed and fitted 
#' frequencies respectively. Both of these statistics converge to 
#' \mjseqn{\chi^2} distribution asymptotically with the same 
#' degrees of freedom.
#' 
#' The convergence of \mjseqn{G, \chi^2} statistics to 
#' \mjseqn{\chi^2} distribution may be violated if expected counts 
#' in cells are too low, say < 5, so it is customary to either censor or 
#' omit these cells.
#' 
#' @param object object of singleRmargin class.
#' @param df degrees of freedom if not provided the function will try and manually
#' but it is not always possible.
#' @param dropl5 a character indicating treatment of cells with frequencies < 5 
#' either grouping them, dropping or leaving them as is. Defaults to drop.
#' @param ... currently does nothing.
#'
#' @method summary singleRmargin
#' @return A chi squared test and G test for comparison between fitted and 
#' observed marginal frequencies.
#' @examples 
#' # Create a simple model
#' Model <- estimatePopsize(
#'   formula = capture ~ ., 
#'   data = netherlandsimmigrant, 
#'   model = ztpoisson, 
#'   method = "IRLS"
#' )
#' plot(Model, "rootogram")
#' # We see a considerable lack of fit
#' summary(marginalFreq(Model), df = 1, dropl5 = "group")
#' @exportS3Method
summary.singleRmargin <- function(object, df,
                                  dropl5 = c("drop", 
                                             "group", 
                                             "no"), 
                                  ...) {
  if (!is.character(dropl5) | length(dropl5) > 1) {
    warning("The argument dropl5 should be a 1 length character vector")
    dropl5 <- dropl5[1]
  }
  
  y <- object$y
  
  if (missing(dropl5) | !is.character(dropl5)) {dropl5 <- "no"}
  if (grepl("zot", object$name) & (1 %in% names(y))) {y <- y[-1]}
  
  if ((missing(df)) && (object$df < 1)) {
    warning("Degrees of freedom may be inacurate.")
    df <- 1
  } else if (missing(df)) {
    df <- object$df
  }
  
  A <- object$table[names(y)]
  
  switch (dropl5,
    "group" = {
      l <- (A < 5)
      if(!all(l == FALSE)) {
        y <- c(y[!l], sum(y[l]))
        A <- c(A[!l], sum(A[l]))
      }
    },
    "drop" = {
      l <- (A < 5)
      y <- y[!l]
      A <- A[!l]
    }
  )
  
  X2 <- sum(((A - y) ^ 2) / A)
  G <- 2 * sum(y * log(y / A))
  
  pval <- stats::pchisq(q = c(X2, G), df = df, lower.tail = FALSE)
  
  vect <- data.frame(
    round(c(X2, G), digits = 2),
    rep(df, 2), signif(pval, digits = 2)
  )
  rownames(vect) <- c("Chi-squared test", "G-test")
  colnames(vect) <- c("Test statistics", "df", "P(>X^2)")
  
  structure(
    list(Test = vect,
         l5 = switch(dropl5,
           "drop"  = "dropped",
           "group" = "grouped",
           "no"    = "preserved"
         ),
         y = y),
    class = "summarysingleRmargin"
  )
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
  invisible()
}