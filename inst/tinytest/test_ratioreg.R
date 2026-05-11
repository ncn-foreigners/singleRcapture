
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

# input validation -----------------------------------------------------------

expect_error(
  ratioReg(
    observed ~ 1,
    data = data.frame(observed = c(1, -1, 2)),
    confint = "none"
  ),
  pattern = "strictly positive"
)

expect_error(
  ratioReg(
    observed ~ 1,
    data = data.frame(observed = c(1, 2.5, 3)),
    confint = "none"
  ),
  pattern = "whole numbers"
)

expect_error(
  ratioReg(
    observed ~ 1,
    data = data.frame(observed = c(1, 2, Inf)),
    confint = "none"
  ),
  pattern = "finite"
)

expect_error(
  ratioReg(
    ~ 1,
    data = data.frame(observed = c(1, 2, 3)),
    confint = "none"
  ),
  pattern = "response variable"
)

expect_error(
  ratioReg(
    observed ~ 0,
    data = data.frame(observed = c(1, 2, 3)),
    confint = "none"
  ),
  pattern = "~ 1"
)

expect_error(
  ratioReg(
    observed ~ 1,
    data = ratio_raw,
    confint = "bootstrap",
    B = -1
  ),
  pattern = "non-negative"
)

expect_error(
  ratioReg(
    observed ~ 1,
    data = ratio_raw,
    confint = "bootstrap",
    B = 0
  ),
  pattern = "at least 1"
)

expect_error(
  ratioReg(
    observed ~ 1,
    data = data.frame(observed = c(rep(1, 10), rep(2, 5))),
    confint = "none",
    maxCount = 1
  ),
  pattern = "at least 2"
)

expect_error(
  ratioReg(
    observed ~ 1,
    data = ratio_raw,
    weights = c(1, -1, rep(1, nrow(ratio_raw) - 2)),
    confint = "none"
  ),
  pattern = "non-negative"
)

# maxCount truncates the support ---------------------------------------------

mc_counts <- data.frame(
  observed = c(rep(1, 100), rep(2, 80), rep(3, 48), rep(4, 19), rep(5, 5))
)
mc_full <- ratioReg(
  observed ~ 1,
  data = mc_counts,
  model = "M0",
  confint = "none"
)
mc_capped <- ratioReg(
  observed ~ 1,
  data = mc_counts,
  model = "M0",
  confint = "none",
  maxCount = 3
)

expect_identical(
  names(mc_full$counts),
  as.character(1:5)
)
expect_identical(
  names(mc_capped$counts),
  as.character(1:3)
)
expect_identical(mc_capped$ratioData$k, c(1L, 2L))

# subset works with both a bare expression and a precomputed logical ---------

set.seed(11)
subset_data <- data.frame(
  observed = sample(1:5, 300, replace = TRUE, prob = c(.4, .3, .15, .1, .05)),
  grp = sample(c("a", "b"), 300, replace = TRUE)
)
keep_a <- subset_data$grp == "a"

fit_bare <- ratioReg(
  observed ~ 1,
  data = subset_data,
  model = "M0",
  confint = "none",
  subset = grp == "a"
)
fit_pre <- ratioReg(
  observed ~ 1,
  data = subset_data,
  model = "M0",
  confint = "none",
  subset = keep_a
)

expect_equal(fit_bare$sizeObserved, sum(keep_a))
expect_equal(
  unname(coef(fit_bare)),
  unname(coef(fit_pre)),
  tolerance = sqrt(.Machine$double.eps)
)

# method variants -----------------------------------------------------------

method_fit <- ratioReg(
  observed ~ 1,
  data = ratio_raw,
  model = "auto",
  estimator = "SM",
  confint = "none"
)

expect_equal(
  nobs(method_fit),
  nrow(method_fit$ratioData)
)

expect_equal(
  unname(extractAIC(method_fit)),
  unname(AIC(method_fit)),
  tolerance = sqrt(.Machine$double.eps)
)

expect_equal(
  unname(coef(method_fit, model = "M0")),
  unname(coef(method_fit, model = "M0", component = "full")),
  tolerance = sqrt(.Machine$double.eps)
)

m1_base <- coef(method_fit, model = "M1", component = "base")
expect_false("oneInflation" %in% names(m1_base))

expect_equal(
  popSizeEst(method_fit)$pointEstimate,
  popSizeEst(method_fit, estimator = method_fit$primaryEstimator)$pointEstimate,
  tolerance = sqrt(.Machine$double.eps)
)

expect_true(
  popSizeEst(method_fit, estimator = "HT")$pointEstimate >=
    method_fit$sizeObserved
)
expect_true(
  popSizeEst(method_fit, estimator = "SM")$pointEstimate >=
    method_fit$sizeObserved
)

# fitted frequencies sum to the observed total ------------------------------

expect_equal(
  unname(sum(method_fit$fittedFrequencies)),
  method_fit$sizeObserved,
  tolerance = 1e-8
)

# fitted base probabilities form a proper distribution ----------------------

expect_equal(
  unname(sum(method_fit$fittedProbabilities)),
  1,
  tolerance = 1e-10
)
expect_true(all(method_fit$fittedProbabilities >= 0))

# custom ratioFormula -------------------------------------------------------

custom_fit <- ratioReg(
  observed ~ 1,
  data = ratio_raw,
  ratioFormula = ~ k,
  model = "M0",
  confint = "none"
)
expect_identical(
  names(coef(custom_fit, model = "M0")),
  c("(Intercept)", "k")
)

# bootstrap with model = "M0" exercises the f0-imputation branch ------------

set.seed(0)
boot_m0 <- ratioReg(
  observed ~ 1,
  data = ratio_raw,
  model = "M0",
  estimator = "HT",
  confint = "bootstrap",
  B = 15,
  seed = 17
)
expect_true(length(popSizeEst(boot_m0)$boot) == 15)
expect_true(is.finite(popSizeEst(boot_m0)$variance))
