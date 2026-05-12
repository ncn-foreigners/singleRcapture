source("helper_families.R")
if (!developer_tests_enabled()) exit_file("developer tests disabled")

inflated_test_data <- load_inflated_test_data()
N <- inflated_test_data$N
x1 <- inflated_test_data$x1
test_inflated <- inflated_test_data$test_inflated

# ztHurdle ####

# poisson
coefficients_ztHurdlepoisson <- c(.6, -.3, -.5, .1)
#eta <- cbind(coefficients_ztHurdlepoisson[1] + coefficients_ztHurdlepoisson[2] * x1, coefficients_ztHurdlepoisson[3] + coefficients_ztHurdlepoisson[4] * x1)

model_ztHurdlepoisson <- ztHurdlepoisson(
  lambdaLink = "log",
  piLink  = "logit",
)

#y <- model_ztHurdlepoisson$simulate(n = N, eta = eta, lower = -1)
observed_ztHurdlepoisson <- test_inflated$ztHurdlepoisson1

data_ztHurdlepoisson <- data.frame(
  x = x1[observed_ztHurdlepoisson > 0],
  y = observed_ztHurdlepoisson[observed_ztHurdlepoisson > 0]
)

fit_ztHurdlepoisson <- estimatePopsize(
  formula = y ~ x,
  model   = model_ztHurdlepoisson,
  data    = data_ztHurdlepoisson,
  method = "IRLS",
  controlModel = controlModel(
    piFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_ztHurdlepoisson))

expect_silent(print(summary(fit_ztHurdlepoisson)))

expect_equivalent(
  coef(fit_ztHurdlepoisson),
  coefficients_ztHurdlepoisson,
  tolerance = .7
)

expect_silent(
  population_ztHurdlepoisson <- popSizeEst(fit_ztHurdlepoisson)
)

expect_equivalent(
  population_ztHurdlepoisson$pointEstimate,
  N,
  tolerance = .5
)

expect_true(
  (population_ztHurdlepoisson$confidenceInterval[1, 1] < N) &
  (N < population_ztHurdlepoisson$confidenceInterval[1, 2])
)

expect_true(
  (population_ztHurdlepoisson$confidenceInterval[2, 1] < N) &
  (N < population_ztHurdlepoisson$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_ztHurdlepoisson,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztHurdlepoisson,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztHurdlepoisson,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_ztHurdlepoisson,
    type = "contr"
  )),
  fit_ztHurdlepoisson$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztHurdlepoisson,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_ztHurdlepoisson,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztHurdlepoisson,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztHurdlepoisson,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztHurdlepoisson,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_ztHurdlepoisson,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)

# same with different link
model_ztHurdlepoisson <- ztHurdlepoisson(
  lambdaLink = "neglog",
  piLink  = "probit",
)

#y <- model_ztHurdlepoisson$simulate(n = N, eta = eta, lower = -1)
observed_ztHurdlepoisson <- test_inflated$ztHurdlepoisson2

data_ztHurdlepoisson <- data.frame(
  x = x1[observed_ztHurdlepoisson > 0],
  y = observed_ztHurdlepoisson[observed_ztHurdlepoisson > 0]
)

fit_ztHurdlepoisson <- estimatePopsize(
  formula = y ~ x,
  model   = model_ztHurdlepoisson,
  data    = data_ztHurdlepoisson,
  method = "IRLS",
  controlModel = controlModel(
    piFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_ztHurdlepoisson))

expect_silent(print(summary(fit_ztHurdlepoisson)))

expect_equivalent(
  coef(fit_ztHurdlepoisson),
  coefficients_ztHurdlepoisson,
  tolerance = .4
)

expect_silent(
  population_ztHurdlepoisson <- popSizeEst(fit_ztHurdlepoisson)
)

expect_equivalent(
  population_ztHurdlepoisson$pointEstimate,
  N,
  tolerance = .4
)

expect_true(
  (population_ztHurdlepoisson$confidenceInterval[1, 1] < N) &
  (N < population_ztHurdlepoisson$confidenceInterval[1, 2])
)

expect_true(
  (population_ztHurdlepoisson$confidenceInterval[2, 1] < N) &
  (N < population_ztHurdlepoisson$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_ztHurdlepoisson,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztHurdlepoisson,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztHurdlepoisson,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_ztHurdlepoisson,
    type = "contr"
  )),
  fit_ztHurdlepoisson$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztHurdlepoisson,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_ztHurdlepoisson,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztHurdlepoisson,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztHurdlepoisson,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztHurdlepoisson,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_ztHurdlepoisson,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)
