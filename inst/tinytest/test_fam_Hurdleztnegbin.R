source("helper_families.R")
if (!developer_tests_enabled()) exit_file("developer tests disabled")

inflated_test_data <- load_inflated_test_data()
N <- inflated_test_data$N
x1 <- inflated_test_data$x1
test_inflated <- inflated_test_data$test_inflated

# negbin
coefficients_Hurdleztnegbin <- c(.4, .4, -.4, .2, -1.1, .37)
#eta <- cbind(coefficients_Hurdleztnegbin[1] + coefficients_Hurdleztnegbin[2] * x1, coefficients_Hurdleztnegbin[3] + coefficients_Hurdleztnegbin[4] * x1, coefficients_Hurdleztnegbin[5] + coefficients_Hurdleztnegbin[6] * x1)

model_Hurdleztnegbin <- Hurdleztnegbin(
  lambdaLink = "log",
  alphaLink  = "log",
  piLink  = "logit",
  eimStep = 30
)

#y <- model_Hurdleztnegbin$simulate(n = N, eta = eta, lower = -1)
observed_Hurdleztnegbin <- test_inflated$Hurdleztnegbin1

data_Hurdleztnegbin <- data.frame(
  x = x1[observed_Hurdleztnegbin > 0],
  y = observed_Hurdleztnegbin[observed_Hurdleztnegbin > 0]
)

fit_Hurdleztnegbin <- estimatePopsize(
  formula = y ~ x,
  model   = model_Hurdleztnegbin,
  data    = data_Hurdleztnegbin,
  method = "IRLS",
  controlModel = controlModel(
    piFormula = ~ x,
    alphaFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_Hurdleztnegbin, confint = TRUE, correlation = TRUE))

expect_silent(print(summary(fit_Hurdleztnegbin, confint = TRUE, correlation = TRUE)))

expect_equivalent(
  coef(fit_Hurdleztnegbin),
  coefficients_Hurdleztnegbin,
  tolerance = .7
)

expect_silent(
  population_Hurdleztnegbin <- popSizeEst(fit_Hurdleztnegbin)
)

expect_equivalent(
  population_Hurdleztnegbin$pointEstimate,
  N,
  tolerance = .35
)

expect_true(
  (population_Hurdleztnegbin$confidenceInterval[2, 1] < N) &
  (N < population_Hurdleztnegbin$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_Hurdleztnegbin,
    type = "response",
    se.fit = TRUE
  )[, 4:6] > 0)
)

expect_true(
  all(predict(
    fit_Hurdleztnegbin,
    type = "link",
    se.fit = TRUE
  )[, 4:6]> 0)
)

expect_true(
  all(predict(
    fit_Hurdleztnegbin,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_Hurdleztnegbin,
    type = "contr"
  )),
  fit_Hurdleztnegbin$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_Hurdleztnegbin,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_Hurdleztnegbin,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 4:6] > 0)
)

expect_true(
  all(predict(
    fit_Hurdleztnegbin,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[, 4:6]> 0)
)

expect_true(
  all(predict(
    fit_Hurdleztnegbin,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_Hurdleztnegbin,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_Hurdleztnegbin,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)
