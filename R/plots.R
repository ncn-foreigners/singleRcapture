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
#' @importFrom graphics boxplot
#' @importFrom graphics panel.smooth
#' @export
plot.singleR <- function(x, 
                         plotType = c("QQ", "marginal", "fitresid",
                                      "bootHist", "rootogram", "dfpopContr",
                                      "dfpopBox", "scaleLoc", "Cooks"), 
                         ...) {
  if ((plotType == "bootHist") && (!is.numeric(x$populationSize$boot))) {
    stop("Trying to plot bootstrap results with no bootstrap performed")
  } 
  plotType <- match.arg(plotType)
  # for now only pearson
  type <- "pearson"
  res <- residuals.singleR(x, type = type)[, 1]
  switch(plotType,
  QQ = {
    stats::qqnorm(res);
    stats::qqline(res);
  },
  marginal = {
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
    FF <- M$table;
    FF[names(M$y)] <- M$y;
    FF[setdiff(names(M$table), names(M$y))] <- 0;
    FF <- FF[-1];
    bp <- graphics::barplot(
      sqrt(FF),
      offset = sqrt(M$table[-1]) - sqrt(FF),
      ylab = expression(sqrt("Frequency")), # This looks just ever so slightly fancier 
      xlab = "captures",
      ylim = c(min(sqrt(M$table[-1]) - sqrt(FF)) - 1, max(sqrt(M$table[-1]) + 1)),
      ...);
    graphics::lines(bp, sqrt(M$table[-1]), 
          type = "o", 
          pch = 19, 
          lwd = 2, 
          col = 2);
    graphics::abline(h = 0, 
           lty = 2)
  },
  "dfpopContr" = {
    dfpop <- dfpopsize(x, observedPop = FALSE, ...);
    contr <- x$model$pointEst(pw = x$prior.weights, 
                              disp = x$dispersion, 
                              lambda = if (x$model$family == "zelterman") {x$model$linkinv(x$X %*% x$coefficients)} else if (grepl(x = x$model$family, pattern = "negbin")) {x$model$linkinv(x$X[rownames(x$linear.predictors), ] %*% x$coefficients[-1])} else {x$model$linkinv(x$X[rownames(x$linear.predictors), ] %*% x$coefficients)}, 
                              contr = TRUE);
    plot(x = dfpop, y = contr,
         main = "Observation deletion effect on point estimate of\npopulation size estimate vs observation contribution",
         xlab = "Deletion effect", ylab = "Observation contribution", 
         ...);
    abline(a = 0, b = 1, col = "red")
  },
  "dfpopBox" = {
    dfpop <- dfpopsize(x, ...);
    graphics::boxplot(dfpop, ylab = "Deletion effect",
                      main = "Boxplot of observation deletion effect on\npoint estimate of population size estimate", 
                      ...)
  },
  "scaleLoc" = {
    plot(y = sqrt(abs(res)), x = x$linear.predictors,
         xlab = "Linear predictors",
         ylab = expression(sqrt("Std. Pearson resid.")),
         main = "Scale-Location plot",
         ...);
    graphics::panel.smooth(y = sqrt(abs(res)), x = x$linear.predictors, iter = 0)
  },
  "fitresid" = {
    plot(y = res, x = x$linear.predictors,
         xlab = "Linear predictors",
         ylab = "Std. Pearson resid.",
         main = "Residuals vs Fitted",
         ...);
    abline(lty = 2, col = "darkgrey", h = 0);
    graphics::panel.smooth(y = res, x = x$linear.predictors, iter = 0)
  },
  "Cooks" = {
    A <- cooks.distance.singleR(x);
    plot(A,
         main = "Cook's distance",
         ylab = "Cook's distance",
         xlab = "Observation index",
         ylim = c(0, max(A) * 1.1),
         ...)
  })
}