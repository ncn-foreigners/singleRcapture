
ratio_raw <- data.frame(
  observed = c(rep(1, 50), rep(2, 30), rep(3, 15), rep(4, 6), rep(5, 2))
)

ratio_tab <- aggregate(
  list(weight = rep(1, nrow(ratio_raw))),
  by = list(observed = ratio_raw$observed),
  FUN = sum
)

expect_error(
  ratioReg(
    observed ~ grp,
    data = data.frame(observed = c(1, 2, 2, 3), grp = c("a", "a", "b", "b")),
    confint = "none"
  ),
  pattern = "ratioFormula"
)

expect_silent(
  ratio_fit_raw <- ratioReg(
    observed ~ 1,
    data = ratio_raw,
    model = "M0",
    estimator = "SM",
    confint = "none"
  )
)

expect_silent(
  ratio_fit_tab <- ratioReg(
    observed ~ 1,
    data = ratio_tab,
    weights = ratio_tab$weight,
    model = "M0",
    estimator = "SM",
    confint = "none"
  )
)

expect_equal(
  coef(ratio_fit_raw),
  coef(ratio_fit_tab),
  tolerance = sqrt(.Machine$double.eps)
)

expect_equal(
  popSizeEst(ratio_fit_raw)$pointEstimate,
  popSizeEst(ratio_fit_tab)$pointEstimate,
  tolerance = sqrt(.Machine$double.eps)
)

gap_counts <- data.frame(
  observed = c(rep(1, 100), rep(2, 80), rep(4, 20), rep(5, 4))
)

gap_fit <- ratioReg(
  observed ~ 1,
  data = gap_counts,
  model = "M0",
  estimator = "HT",
  confint = "none"
)

expect_identical(
  gap_fit$ratioData$k,
  c(1L, 4L)
)

contig_counts <- data.frame(
  observed = c(rep(1, 100), rep(2, 80), rep(3, 48), rep(4, 19))
)

contig_fit <- ratioReg(
  observed ~ 1,
  data = contig_counts,
  model = "M0",
  estimator = "HT",
  confint = "none"
)

manual_ratio_data <- data.frame(
  k = 1:3,
  ratio = c(80 / 100, 48 / 80, 19 / 48)
)
manual_ratio_data$logRatio <- log(manual_ratio_data$ratio)
manual_ratio_data$weight <- 1 / c(
  1 / 100 + 1 / 80,
  1 / 80 + 1 / 48,
  1 / 48 + 1 / 19
)

manual_lm <- lm(
  logRatio ~ log(k + 1),
  data = manual_ratio_data,
  weights = weight
)

expect_equal(
  coef(contig_fit),
  coef(manual_lm),
  tolerance = 1e-10
)

manual_sse <- sum(manual_ratio_data$weight * residuals(manual_lm) ^ 2)
manual_sigma2 <- manual_sse / nrow(manual_ratio_data)
manual_loglik <- -0.5 * nrow(manual_ratio_data) *
  (log(2 * pi) + 1 + log(manual_sigma2)) +
  0.5 * sum(log(manual_ratio_data$weight))

expect_equal(
  as.numeric(logLik(contig_fit)),
  manual_loglik,
  tolerance = 1e-10
)

expect_equal(
  AIC(contig_fit),
  -2 * manual_loglik + 2 * (manual_lm$rank + 1),
  tolerance = 1e-10
)

expect_equal(
  BIC(contig_fit),
  -2 * manual_loglik + log(nrow(manual_ratio_data)) * (manual_lm$rank + 1),
  tolerance = 1e-10
)

manual_log_ratios <- coef(contig_fit)[1] + coef(contig_fit)[2] * log(1:4)
manual_terms <- exp(c(0, cumsum(manual_log_ratios)))
manual_p0 <- 1 / sum(manual_terms)
manual_ht <- sum(contig_fit$counts) / (1 - manual_p0)
manual_sm <- sum(contig_fit$counts) + contig_fit$counts["2"] /
  exp(sum(manual_log_ratios[1:2]))

expect_equal(
  unname(contig_fit$fittedProbabilities["0"]),
  manual_p0,
  tolerance = 1e-10
)

expect_equal(
  popSizeEst(contig_fit, estimator = "HT")$pointEstimate,
  manual_ht,
  tolerance = 1e-10
)

expect_equal(
  unname(popSizeEst(contig_fit, estimator = "SM")$pointEstimate),
  unname(manual_sm),
  tolerance = 1e-10
)

no_infl_counts <- data.frame(
  observed = c(rep(1, 30), rep(2, 27), rep(3, 20), rep(4, 12), rep(5, 7), rep(6, 4))
)

no_infl_fit <- ratioReg(
  observed ~ 1,
  data = no_infl_counts,
  model = "auto",
  estimator = "SM",
  confint = "none"
)

inflated_counts <- data.frame(
  observed = c(rep(1, 200), rep(2, 30), rep(3, 15), rep(4, 6), rep(5, 2))
)

inflated_fit <- ratioReg(
  observed ~ 1,
  data = inflated_counts,
  model = "auto",
  estimator = "SM",
  confint = "none"
)

expect_identical(
  no_infl_fit$selectedModel,
  "M0"
)

expect_identical(
  inflated_fit$selectedModel,
  "M1"
)

expect_silent(
  boot_fit <- ratioReg(
    observed ~ 1,
    data = ratio_raw,
    model = "auto",
    estimator = "SM",
    confint = "bootstrap",
    B = 20,
    seed = 123
  )
)

expect_silent(
  boot_fit_same <- ratioReg(
    observed ~ 1,
    data = ratio_raw,
    model = "auto",
    estimator = "SM",
    confint = "bootstrap",
    B = 20,
    seed = 123
  )
)

expect_equal(
  popSizeEst(boot_fit)$confidenceInterval,
  popSizeEst(boot_fit_same)$confidenceInterval,
  tolerance = sqrt(.Machine$double.eps)
)

expect_equal(
  confint(boot_fit),
  as.matrix(popSizeEst(boot_fit)$confidenceInterval),
  tolerance = sqrt(.Machine$double.eps)
)

expect_true(
  length(popSizeEst(boot_fit)$boot) == 20
)

expect_silent(summary(boot_fit))
expect_silent(print(boot_fit))

expect_silent({
  grDevices::pdf(tempfile(fileext = ".pdf"))
  plot(boot_fit, plotType = "ratio")
  grDevices::dev.off()
})

expect_silent({
  grDevices::pdf(tempfile(fileext = ".pdf"))
  plot(boot_fit, plotType = "bootHist")
  grDevices::dev.off()
})

expect_silent({
  grDevices::pdf(tempfile(fileext = ".pdf"))
  plot(boot_fit, plotType = "freqFit")
  grDevices::dev.off()
})

expect_error(
  estimatePopsize(
    capture ~ 1,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    ratioReg = TRUE,
    controlMethod = controlMethod(silent = TRUE)
  ),
  pattern = "Use ratioReg"
)
