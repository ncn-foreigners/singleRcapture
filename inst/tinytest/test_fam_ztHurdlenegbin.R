source("helper_families.R")
if (!developer_tests_enabled()) exit_file("developer tests disabled")

inflated_test_data <- load_inflated_test_data()
N <- inflated_test_data$N
x1 <- inflated_test_data$x1
test_inflated <- inflated_test_data$test_inflated

# negbin
coefficients_ztHurdlenegbin <- c(.4, .4, -.4, .2, -1.1, .37)
#eta <- cbind(coefficients_ztHurdlenegbin[1] + coefficients_ztHurdlenegbin[2] * x1, coefficients_ztHurdlenegbin[3] + coefficients_ztHurdlenegbin[4] * x1, coefficients_ztHurdlenegbin[5] + coefficients_ztHurdlenegbin[6] * x1)

model_ztHurdlenegbin <- ztHurdlenegbin(
  lambdaLink = "log",
  alphaLink  = "log",
  piLink  = "logit",
  eimStep = 20
)

#y <- model_ztHurdlenegbin$simulate(n = N, eta = eta, lower = -1)
observed_ztHurdlenegbin <- test_inflated$ztHurdlenegbin1

data_ztHurdlenegbin <- data.frame(
  x = x1[observed_ztHurdlenegbin > 0],
  y = observed_ztHurdlenegbin[observed_ztHurdlenegbin > 0]
)

fit_ztHurdlenegbin <- estimatePopsize(
  formula = y ~ x,
  model   = model_ztHurdlenegbin,
  data    = data_ztHurdlenegbin,
  method = "IRLS",
  controlModel = controlModel(
    piFormula = ~ x,
    alphaFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_ztHurdlenegbin, confint = TRUE, correlation = TRUE))

expect_silent(print(summary(fit_ztHurdlenegbin, confint = TRUE, correlation = TRUE)))

expect_equivalent(
  coef(fit_ztHurdlenegbin),
  coefficients_ztHurdlenegbin,
  tolerance = .5
)

expect_silent(
  population_ztHurdlenegbin <- popSizeEst(fit_ztHurdlenegbin)
)

expect_equivalent(
  population_ztHurdlenegbin$pointEstimate,
  N,
  tolerance = .3
)

expect_true(
  (population_ztHurdlenegbin$confidenceInterval[1, 1] < N) &
  (N < population_ztHurdlenegbin$confidenceInterval[1, 2])
)

expect_true(
  (population_ztHurdlenegbin$confidenceInterval[2, 1] < N) &
  (N < population_ztHurdlenegbin$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_ztHurdlenegbin,
    type = "response",
    se.fit = TRUE
  )[, 4:6] > 0)
)

expect_true(
  all(predict(
    fit_ztHurdlenegbin,
    type = "link",
    se.fit = TRUE
  )[, 4:6]> 0)
)

expect_true(
  all(predict(
    fit_ztHurdlenegbin,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_ztHurdlenegbin,
    type = "contr"
  )),
  fit_ztHurdlenegbin$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztHurdlenegbin,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_ztHurdlenegbin,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 4:6] > 0)
)

expect_true(
  all(predict(
    fit_ztHurdlenegbin,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[, 4:6]> 0)
)

expect_true(
  all(predict(
    fit_ztHurdlenegbin,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztHurdlenegbin,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_ztHurdlenegbin,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)
