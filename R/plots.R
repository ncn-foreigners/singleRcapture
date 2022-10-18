#' @title Diagnostic plots singleR regression
#'
#' @param x object of singleRmargin class.
#' @param plotType Character parameter specifying type of plot to be made.
#' The following list presents and briefly explains possible type of plots:
#' \itemize{
#'   \item \code{QQ} --
#'   \item \code{marginal} --
#'   \item \code{fitresid} --
#'   \item \code{bootHist} --
#'   \item \code{rootogram} --
#'   \item \code{dfpopContr} --
#'   \item \code{dfpopBox} --
#'   \item \code{scaleLoc} --
#'   \item \code{Cooks} --
#'   \item \code{hatplot} --
#' }
#' @param ... 
#' \itemize{
#'   \item 
#'   \item
#'   \item
#'   \item
#' }
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
                                      "dfpopBox", "scaleLoc", "Cooks",
                                      "hatplot"), 
                         ...) {
  if ((plotType == "bootHist") && (!is.numeric(x$populationSize$boot))) {
    stop("Trying to plot bootstrap results with no bootstrap performed")
  } 
  plotType <- match.arg(plotType)
  # move this to particular plots
  if (x$model$parNum == 1) type <- "pearsonSTD" else type <- "pearson"
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
  dfpopContr = {
    dfpop <- dfpopsize(x, observedPop = if (x$model$family == "zelterman") TRUE else FALSE, ...);
    contr <- x$model$pointEst(pw = x$prior.weights[x$which$est], 
                              eta = x$linear.predictors,
                              contr = TRUE);
    plot(x = dfpop, y = contr,
         main = "Observation deletion effect on point estimate of\npopulation size estimate vs observation contribution",
         xlab = "Deletion effect", ylab = "Observation contribution", 
         ...);
    abline(a = 0, b = 1, col = "red")
  },
  dfpopBox = {
    dfpop <- dfpopsize(x, observedPop = FALSE,...);
    graphics::boxplot(dfpop, ylab = "Deletion effect",
                      main = "Boxplot of observation deletion effect on\npoint estimate of population size estimate", 
                      ...)
  },
  scaleLoc = {
    if (x$model$parNum == 1) {
      plot(y = sqrt(abs(res)), x = x$linear.predictors,
           xlab = "Linear predictors",
           ylab = expression(sqrt("Std. Pearson resid.")),
           main = "Scale-Location plot",
           ...);
      graphics::panel.smooth(y = sqrt(abs(res)), x = x$linear.predictors, iter = 0)
    } else {
      par <- graphics::par(no.readonly = TRUE);
      par(mfrow = c(x$model$parNum, 1));
      for (k in 1:x$model$parNum) {
        plot(y = sqrt(abs(res)), x = x$linear.predictors[, k],
             xlab = "Linear predictors",
             ylab = expression(sqrt("Pearson resid.")),
             main = "Scale-Location plot",
             sub = paste0("For linear predictors associated with: ", x$model$etaNames[k]),
             ...);
        graphics::panel.smooth(y = sqrt(abs(res)), x = x$linear.predictors[, k], iter = 0)
      }
      graphics::par(par);
    }
  },
  fitresid = {
    if (x$model$parNum == 1) {
      plot(y = res, x = x$linear.predictors,
           xlab = "Linear predictors",
           ylab = "Std. Pearson resid.",
           main = "Residuals vs Fitted",
           ...);
      abline(lty = 2, col = "darkgrey", h = 0);
      graphics::panel.smooth(y = res, x = x$linear.predictors, iter = 0)
    } else {
      par <- graphics::par(no.readonly = TRUE);
      par(mfrow = c(x$model$parNum, 1));
      for (k in 1:x$model$parNum) {
        plot(y = res, x = x$linear.predictors[, k],
             xlab = "Linear predictors",
             ylab = "Pearson resid.",
             main = "Residuals vs Fitted",
             sub = paste0("For linear predictors associated with: ", x$model$etaNames[k]),
             ...);
        abline(lty = 2, col = "darkgrey", h = 0);
        graphics::panel.smooth(y = res, x = x$linear.predictors[, k], iter = 0)
      }
      graphics::par(par);
    }
  },
  Cooks = {
    A <- cooks.distance.singleR(x);
    plot(A,
         main = "Cook's distance",
         ylab = "Cook's distance",
         xlab = "Observation index",
         ylim = c(0, max(A) * 1.1),
         ...)
  },
  hatplot = {
    A <- hatvalues.singleR(x, ...);
    par <- graphics::par(no.readonly = TRUE);
    par(mfrow = c(x$model$parNum, 1));
    for (k in 1:x$model$parNum) {
      plot(A[, k],
           xlab = "Observation index",
           ylab = "Hat values",
           main = paste0("For linear predictors associated with: ", x$model$etaNames[k]),
           ...)
    }
    graphics::par(par);
  })
}