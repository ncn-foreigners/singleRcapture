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
#' @importFrom graphics matplot 
#' @importFrom graphics legend
#' @export
plot.singleR <- function(x, 
                         plotType = c("QQ", "marginal",
                                      "fitresid"), 
                         ...) {
  # This is currently mostly a placeholder to remind myself that it needs to be done at some point
  # TODO:
  # QQ - standardardise
  # add all other from plot.lm
  # marginal - optional y, optional 0, optional more models
  plotType <- match.arg(plotType)
  if(missing(type)) {type <- "pearson"}
  res <- residuals.singleR(x, type = type)
  switch(plotType,
  QQ = plot(
    x = sort(res),
    y = stats::qnorm(stats::ppoints(length(res))), 
    abline(a = 0, b = 1 , lty = 2),
    ...),
  marginal = graphics::matplot(
    y = cbind(marginalFreq(x)$table, c("0" = 0, table(x$y))),
    x = 0:max(x$y), type = "o",
    col = 1:2, pch = 21:22, lty = 2,
    main = "Plot of observed and fitted marginal frequencies",
    ylab = "Frequency",
    xlab = "Counts",
    ...))
  if (type == "marginal") {
    legend("topright",
           legend = c(x$model$family,
                      "Observed"),
           col = 1:2, pch = 21:22)
  }
}