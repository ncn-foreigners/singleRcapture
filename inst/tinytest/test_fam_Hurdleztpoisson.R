source("helper_families.R")
if (!developer_tests_enabled()) exit_file("developer tests disabled")

inflated_test_data <- load_inflated_test_data()
N <- inflated_test_data$N
x1 <- inflated_test_data$x1
test_inflated <- inflated_test_data$test_inflated

# Hurdlezt ####
# poisson
coefficients_Hurdleztpoisson <- c(.6, -.3, -.5, .1)
#eta <- cbind(coefficients_Hurdleztpoisson[1] + coefficients_Hurdleztpoisson[2] * x1, coefficients_Hurdleztpoisson[3] + coefficients_Hurdleztpoisson[4] * x1)

model_Hurdleztpoisson <- Hurdleztpoisson(
  lambdaLink = "log",
  piLink  = "logit",
)

#y <- model_Hurdleztpoisson$simulate(n = N, eta = eta, lower = -1)
observed_Hurdleztpoisson <- test_inflated$Hurdleztpoisson1

data_Hurdleztpoisson <- data.frame(
  x = x1[observed_Hurdleztpoisson > 0],
  y = observed_Hurdleztpoisson[observed_Hurdleztpoisson > 0]
)

fit_Hurdleztpoisson <- estimatePopsize(
  formula = y ~ x,
  model   = model_Hurdleztpoisson,
  data    = data_Hurdleztpoisson,
  method = "IRLS",
  controlModel = controlModel(
    piFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_Hurdleztpoisson))

expect_silent(print(summary(fit_Hurdleztpoisson)))

expect_equivalent(
  coef(fit_Hurdleztpoisson),
  coefficients_Hurdleztpoisson,
  tolerance = .6
)

expect_silent(
  population_Hurdleztpoisson <- popSizeEst(fit_Hurdleztpoisson)
)

expect_equivalent(
  population_Hurdleztpoisson$pointEstimate,
  N,
  tolerance = .25
)

expect_true(
  (population_Hurdleztpoisson$confidenceInterval[1, 1] < N) &
  (N < population_Hurdleztpoisson$confidenceInterval[1, 2])
)

expect_true(
  (population_Hurdleztpoisson$confidenceInterval[2, 1] < N) &
  (N < population_Hurdleztpoisson$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_Hurdleztpoisson,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_Hurdleztpoisson,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_Hurdleztpoisson,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_Hurdleztpoisson,
    type = "contr"
  )),
  fit_Hurdleztpoisson$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_Hurdleztpoisson,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_Hurdleztpoisson,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_Hurdleztpoisson,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_Hurdleztpoisson,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_Hurdleztpoisson,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_Hurdleztpoisson,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)

# same with different link
model_Hurdleztpoisson <- Hurdleztpoisson(
  lambdaLink = "neglog",
  piLink  = "probit",
)

#y <- model_Hurdleztpoisson$simulate(n = N, eta = eta, lower = -1)
observed_Hurdleztpoisson <- test_inflated$Hurdleztpoisson2

data_Hurdleztpoisson <- data.frame(
  x = x1[observed_Hurdleztpoisson > 0],
  y = observed_Hurdleztpoisson[observed_Hurdleztpoisson > 0]
)

fit_Hurdleztpoisson <- estimatePopsize(
  formula = y ~ x,
  model   = model_Hurdleztpoisson,
  data    = data_Hurdleztpoisson,
  method = "IRLS",
  controlModel = controlModel(
    piFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE, stepsize = .3)
)

expect_silent(summary(fit_Hurdleztpoisson))

expect_silent(print(summary(fit_Hurdleztpoisson)))

expect_equivalent(
  coef(fit_Hurdleztpoisson),
  coefficients_Hurdleztpoisson,
  tolerance = .6
)

expect_silent(
  population_Hurdleztpoisson <- popSizeEst(fit_Hurdleztpoisson)
)

expect_equivalent(
  population_Hurdleztpoisson$pointEstimate,
  N,
  tolerance = .3
)

expect_true(
  (population_Hurdleztpoisson$confidenceInterval[1, 1] < N) &
  (N < population_Hurdleztpoisson$confidenceInterval[1, 2])
)

expect_true(
  (population_Hurdleztpoisson$confidenceInterval[2, 1] < N) &
  (N < population_Hurdleztpoisson$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_Hurdleztpoisson,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_Hurdleztpoisson,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_Hurdleztpoisson,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_Hurdleztpoisson,
    type = "contr"
  )),
  fit_Hurdleztpoisson$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_Hurdleztpoisson,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_Hurdleztpoisson,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_Hurdleztpoisson,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_Hurdleztpoisson,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_Hurdleztpoisson,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_Hurdleztpoisson,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)
