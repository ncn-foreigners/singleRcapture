source("helper_families.R")
if (!developer_tests_enabled()) exit_file("developer tests disabled")

inflated_test_data <- load_inflated_test_data()
N <- inflated_test_data$N
x1 <- inflated_test_data$x1
test_inflated <- inflated_test_data$test_inflated

# negbin
coefficients_oiztnegbin <- c(.4, .4, -.4, .2, -1.1, .37)
#eta <- cbind(coefficients_oiztnegbin[1] + coefficients_oiztnegbin[2] * x1, coefficients_oiztnegbin[3] + coefficients_oiztnegbin[4] * x1, coefficients_oiztnegbin[5] + coefficients_oiztnegbin[6] * x1)

model_oiztnegbin <- oiztnegbin(
  lambdaLink = "log",
  alphaLink  = "log",
  omegaLink  = "logit",
)

#y <- model_oiztnegbin$simulate(n = N, eta = eta, lower = -1)
observed_oiztnegbin <- test_inflated$oiztnegbin1

data_oiztnegbin <- data.frame(
  x = x1[observed_oiztnegbin > 0],
  y = observed_oiztnegbin[observed_oiztnegbin > 0]
)

fit_oiztnegbin <- estimatePopsize(
  formula = y ~ x,
  model   = model_oiztnegbin,
  data    = data_oiztnegbin,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ x,
    alphaFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_oiztnegbin, confint = TRUE, correlation = TRUE))

expect_silent(print(summary(fit_oiztnegbin, confint = TRUE, correlation = TRUE)))

expect_equivalent(
  coef(fit_oiztnegbin),
  coefficients_oiztnegbin,
  tolerance = .6
)

expect_silent(
  population_oiztnegbin <- popSizeEst(fit_oiztnegbin)
)

expect_equivalent(
  population_oiztnegbin$pointEstimate,
  N,
  tolerance = .2
)

expect_true(
  (population_oiztnegbin$confidenceInterval[1, 1] < N) &
  (N < population_oiztnegbin$confidenceInterval[1, 2])
)

expect_true(
  (population_oiztnegbin$confidenceInterval[2, 1] < N) &
  (N < population_oiztnegbin$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_oiztnegbin,
    type = "response",
    se.fit = TRUE
  )[, 4:6] > 0)
)

expect_true(
  all(predict(
    fit_oiztnegbin,
    type = "link",
    se.fit = TRUE
  )[, 4:6]> 0)
)

expect_true(
  all(predict(
    fit_oiztnegbin,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_oiztnegbin,
    type = "contr"
  )),
  fit_oiztnegbin$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_oiztnegbin,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_oiztnegbin,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 4:6] > 0)
)

expect_true(
  all(predict(
    fit_oiztnegbin,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[, 4:6]> 0)
)

expect_true(
  all(predict(
    fit_oiztnegbin,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_oiztnegbin,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_oiztnegbin,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)
