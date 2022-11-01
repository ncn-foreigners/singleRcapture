#' @title Diagnostic plots for regression and population size estimation.
#'
#' @description Simple diagnostic plots for \code{singleR} class objects.
#'
#' @param x Object of \code{singleR} class.
#' @param plotType Character parameter specifying type of plot to be made.
#' The following list presents and briefly explains possible type of plots:
#' \itemize{
#'   \item \code{QQ} -- The quantile-quantile plot for pearson residuals 
#'   (or standardised pearson residuals if these are available for the model) i.e. 
#'   empirical quantiles from residuals are plotted against theoretical quantiles 
#'   from standard distribution.
#'   \item \code{marginal} -- A plot made by \code{matplot} with fitted and 
#'   observed marginal frequencies with labels.
#'   \item \code{fitresid} -- Plot of fitted linear predictors against 
#'   (standardised) pearson residuals.
#'   \item \code{bootHist} -- Simple histogram for statistics obtained from 
#'   bootstraping (if one was performed and the statistics were saved).
#'   \item \code{rootogram} -- Rootogram, for full explanation see: 
#'   Kleiber and Zeileis Visualizing Count Data Regressions Using Rootograms (2016), 
#'   in short it is a \code{barplot} where height is the square root of observed marginal 
#'   frequencies adjusted by difference between square root of observed and fitted marginal 
#'   frequencies connected by line representing fitted marginal frequencies. 
#'   The less of a difference there is between the 0 line and beginning of a bar 
#'   the more accurate fitt was produced by the model.
#'   \item \code{dfpopContr} -- Plot of \code{dfpopsize} against unit contribution.
#'   On the plot is y = x line i.e. what deletion effect would be if removing the
#'   unit from the model didn't effect regression coefficients. The further away
#'   the observation is from this line the more influential it is.
#'   \item \code{dfpopBox} -- Boxplot of \code{dfpopsize} for getting the general 
#'   idea about the distribution of the "influence" of each unit on
#'   population size estimate.
#'   \item \code{scaleLoc} -- The scale - location plot i.e. square root of 
#'   absolute values of (standardised) pearson residuals against linear predictors
#'   for each column of linear predictors.
#'   \item \code{Cooks} -- Plot of cooks distance for detecting influential observations.
#'   \item \code{hatplot} -- Plot of hat values for each linear predictor for detecting influential observations.
#' }
#' @param ... Additional optional arguments passed to the following functions:
#' \itemize{
#'   \item For \code{plotType = "bootHist"}
#'   \itemize{
#'   \item \code{graphics::hist} -- with \code{x, main, xlab, ylab} parameters fixed.
#'   }
#'   \item For \code{plotType = "rootogram"} 
#'   \itemize{
#'   \item \code{graphics::barplot} -- with \code{height, offset, ylab, xlab, ylim} parameters fixed.
#'   \item \code{graphics::lines} -- with \code{x, y, pch, type, lwd, col} parameters fixed.
#'   }
#'   \item For \code{plotType = "dfpopContr"}
#'   \itemize{
#'   \item \code{dfpopsize} -- with \code{model, observedPop} parameters fixed.
#'   \item \code{plot.default} -- with \code{x, y, xlab, main} parameters fixed.
#'   }
#'   \item For \code{plotType = "dfpopBox"}
#'   \itemize{
#'   \item \code{dfpopsize} -- with \code{model, observedPop} parameters fixed.
#'   \item \code{graphics::boxplot} -- with \code{x, ylab, main} parameters fixed.
#'   }
#'   \item For \code{plotType = "scaleLoc"}
#'   \itemize{
#'   \item \code{plot.default} -- with \code{x, y, xlab, ylab, main, sub} parameters fixed.
#'   }
#'   \item For \code{plotType = "fitresid"}
#'   \itemize{
#'   \item \code{plot.default} -- with \code{x, y, xlab, ylab, main, sub} parameters fixed.
#'   }
#'   \item For \code{plotType = "Cooks"}
#'   \itemize{
#'   \item \code{plot.default} -- with \code{x, xlab, ylab, main} parameters fixed.
#'   }
#'   \item For \code{plotType = "hatplot"} 
#'   \itemize{
#'   \item \code{hatvalues.singleR}
#'   \item \code{plot.default} -- with \code{x, xlab, ylab, main} parameters fixed.
#'   }
#' }
#' 
#' @method plot singleR
#' @return No return value only the plot being made.
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
#' @seealso [estimate_popsize()] [dfpopsize()] [marginalFreq()] [stats::plot.lm()] [stats::cooks.distance()] [hatvalues.singleR()]
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
      lp <- x$linear.predictors
      if (x$model$family == "zelterman") {
        lp <- lp[x$which$reg]
      }
      plot(y = sqrt(abs(res)), x = lp,
           xlab = "Linear predictors",
           ylab = expression(sqrt("Std. Pearson resid.")),
           main = "Scale-Location plot",
           ...);
      graphics::panel.smooth(y = sqrt(abs(res)), x = lp, iter = 0)
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
      lp <- x$linear.predictors
      if (x$model$family == "zelterman") {
        lp <- lp[x$which$reg]
      }
      plot(y = res, x = lp,
           xlab = "Linear predictors",
           ylab = "Std. Pearson resid.",
           main = "Residuals vs Fitted",
           ...);
      abline(lty = 2, col = "darkgrey", h = 0);
      graphics::panel.smooth(y = res, x = lp, iter = 0)
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
  invisible()
}