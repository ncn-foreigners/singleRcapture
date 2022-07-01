#' @title Diagnostic plots singleR regression
#'
#' @param x TODO
#' @param plotType TODO
#' @param ... TODO
#' 
#' @method plot singleR
#' @return TODO
#' @importFrom stats ppoints
#' @importFrom graphics abline
#' @importFrom graphics barplot
#' @importFrom graphics hist
#' @importFrom graphics lines
#' @importFrom graphics matplot 
#' @importFrom graphics legend
#' @importFrom stats qqline
#' @importFrom stats qqnorm
#' @export
plot.singleR <- function(x, 
                         plotType = c("QQ", 
                                      "marginal",
                                      "fitresid",
                                      "bootHist",
                                      "rootogram"), 
                         ...) {
  if (plotType == "bootHist" & (!is.numeric(x$populationSize$boot))) {
    stop("Trying to plot bootstrap results with no bootstrap performed")
  } 
  # This is currently mostly a placeholder to remind myself that it needs to be done at some point
  # TODO:
  # QQ - standardardise
  # add all other from plot.lm
  # marginal - optional y, optional 0, optional more models
  plotType <- match.arg(plotType)
  # for now only pearson
  type <- "pearson"
  res <- residuals.singleR(x, type = type)[, 1]
  switch(plotType,
  QQ = {
    stats::qqnorm(res);
    stats::qqline(res);
  },
  marginal = { # TODO add more options
    M <- marginalFreq(x);
    FF <- M$table;
    FF[names(M$y)] <- M$y;
    FF[setdiff(names(M$table), names(M$y))] <- 0;
    graphics::matplot(
      y = cbind(M$table, FF),
      x = 0:max(x$y), type = "o",
      col = 1:2, pch = 21:22, lty = 2,
      main = "Plot of observed and fitted marginal frequencies",
      ylab = "Frequency",
      xlab = "Counts",
      ...);
    legend("topright",
           legend = c(x$model$family,
                      "Observed"),
           col = 1:2, pch = 21:22)
  },
  bootHist = graphics::hist(
    x$populationSize$boot,
    ylab = "Number of bootstrap samples",
    xlab = expression(hat(N)),
    main = "Bootstrap of population size estimates",
    ...
  ),
  rootogram = {
    M <- marginalFreq(x);
    bp <- graphics::barplot(
      sqrt(M$y),
      offset = sqrt(M$table[-1]) - sqrt(M$y),
      ylab = "sqrt(Frequency)", 
      xlab = "captures",
      ...);
    graphics::lines(bp, sqrt(M$table[-1]), 
          type = "o", 
          pch = 19, 
          lwd = 2, 
          col = 2);
    graphics::abline(h = 0, 
           lty = 2)
  }
  )
}