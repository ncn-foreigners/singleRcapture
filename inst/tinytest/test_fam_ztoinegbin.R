source("helper_families.R")
if (!developer_tests_enabled()) exit_file("developer tests disabled")

inflated_test_data <- load_inflated_test_data()
N <- inflated_test_data$N
x1 <- inflated_test_data$x1
test_inflated <- inflated_test_data$test_inflated

# ztoinegbin
coefficients_ztoinegbin <- c(.6, .2, -1.25, .1, -1, -.3)
#eta <- cbind(coefficients_ztoinegbin[1] + coefficients_ztoinegbin[2] * x1, coefficients_ztoinegbin[3] + coefficients_ztoinegbin[4] * x1, coefficients_ztoinegbin[5] + coefficients_ztoinegbin[6] * x1)

model_ztoinegbin <- ztoinegbin(
  lambdaLink = "log",
  alphaLink  = "log",
  omegaLink  = "logit",
)

#y <- model_ztoinegbin$simulate(n = N, eta = eta, lower = -1)
observed_ztoinegbin <- test_inflated$ztoinegbin1

data_ztoinegbin <- data.frame(
  x = x1[observed_ztoinegbin > 0],
  y = observed_ztoinegbin[observed_ztoinegbin > 0]
)

fit_ztoinegbin <- estimatePopsize(
  formula = y ~ x,
  model   = model_ztoinegbin,
  data    = data_ztoinegbin,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ x,
    alphaFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_ztoinegbin))

expect_silent(print(summary(fit_ztoinegbin)))

# don't change this
expect_true(
  sum(abs(coef(fit_ztoinegbin) - coefficients_ztoinegbin)) < 8
)

expect_silent(
  population_ztoinegbin <- popSizeEst(fit_ztoinegbin)
)

expect_equivalent(
  population_ztoinegbin$pointEstimate,
  N,
  tolerance = .3
)

expect_true(
  (population_ztoinegbin$confidenceInterval[1, 1] < N) &
  (N < population_ztoinegbin$confidenceInterval[1, 2])
)

expect_true(
  (population_ztoinegbin$confidenceInterval[2, 1] < N) &
  (N < population_ztoinegbin$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_ztoinegbin,
    type = "response",
    se.fit = TRUE
  )[, 4:6] > 0)
)

expect_true(
  all(predict(
    fit_ztoinegbin,
    type = "link",
    se.fit = TRUE
  )[, 4:6]> 0)
)

expect_true(
  all(predict(
    fit_ztoinegbin,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_ztoinegbin,
    type = "contr"
  )),
  fit_ztoinegbin$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztoinegbin,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_ztoinegbin,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 4:6] > 0)
)

expect_true(
  all(predict(
    fit_ztoinegbin,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[, 4:6]> 0)
)

expect_true(
  all(predict(
    fit_ztoinegbin,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztoinegbin,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_ztoinegbin,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)
