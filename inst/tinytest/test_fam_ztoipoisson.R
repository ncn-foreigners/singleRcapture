source("helper_families.R")
if (!developer_tests_enabled()) exit_file("developer tests disabled")

inflated_test_data <- load_inflated_test_data()
N <- inflated_test_data$N
x1 <- inflated_test_data$x1
test_inflated <- inflated_test_data$test_inflated

# Ztoi ####

# ztoipoisson
coefficients_ztoipoisson <- c(.6, .2, -1.25, .1)
#eta <- cbind(coefficients_ztoipoisson[1] + coefficients_ztoipoisson[2] * x1, coefficients_ztoipoisson[3] + coefficients_ztoipoisson[4] * x1)

model_ztoipoisson <- ztoipoisson(
  lambdaLink = "log",
  omegaLink  = "logit",
)

#y <- model_ztoipoisson$simulate(n = N, eta = eta, lower = -1)
observed_ztoipoisson <- test_inflated$ztoipoisson1

data_ztoipoisson <- data.frame(
  x = x1[observed_ztoipoisson > 0],
  y = observed_ztoipoisson[observed_ztoipoisson > 0]
)

fit_ztoipoisson <- estimatePopsize(
  formula = y ~ x,
  model   = model_ztoipoisson,
  data    = data_ztoipoisson,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_ztoipoisson))

expect_silent(print(summary(fit_ztoipoisson)))

expect_equivalent(
  coef(fit_ztoipoisson),
  coefficients_ztoipoisson, tol = .3
)

expect_silent(
  population_ztoipoisson <- popSizeEst(fit_ztoipoisson)
)

expect_equivalent(
  population_ztoipoisson$pointEstimate,
  N, tol = .1
)

expect_true(
  (population_ztoipoisson$confidenceInterval[1, 1] < N) &
  (N < population_ztoipoisson$confidenceInterval[1, 2])
)

expect_true(
  (population_ztoipoisson$confidenceInterval[2, 1] < N) &
  (N < population_ztoipoisson$confidenceInterval[2, 2])
)

# test whether estimation of subpopulation sizes are 'good'
strata_ztoipoisson <- stratifyPopsize(fit_ztoipoisson, strata = ~ x)
strata_order_ztoipoisson <- 1 + as.numeric(substr(strata_ztoipoisson$name, start = 4, stop = 4))

expect_equivalent(
  strata_ztoipoisson$Estimated,
  as.numeric(table(x1))[strata_order_ztoipoisson],
  tolerance = .15
)

expect_true(
  all(predict(
    fit_ztoipoisson,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztoipoisson,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztoipoisson,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_ztoipoisson,
    type = "contr"
  )),
  fit_ztoipoisson$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztoipoisson,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_ztoipoisson,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztoipoisson,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztoipoisson,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztoipoisson,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_ztoipoisson,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)

# same with different link
model_ztoipoisson <- ztoipoisson(
  lambdaLink = "neglog",
  omegaLink  = "cloglog",
)

#y <- model_ztoipoisson$simulate(n = N, eta = eta, lower = -1)
observed_ztoipoisson <- test_inflated$ztoipoisson2

data_ztoipoisson <- data.frame(
  x = x1[observed_ztoipoisson > 0],
  y = observed_ztoipoisson[observed_ztoipoisson > 0]
)

fit_ztoipoisson <- estimatePopsize(
  formula = y ~ x,
  model   = model_ztoipoisson,
  data    = data_ztoipoisson,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_ztoipoisson))

expect_silent(print(summary(fit_ztoipoisson)))

expect_equivalent(
  coef(fit_ztoipoisson),
  coefficients_ztoipoisson,
  tolerance = .7
)

expect_silent(
  population_ztoipoisson <- popSizeEst(fit_ztoipoisson)
)

# expect_equivalent(
#   population_ztoipoisson$pointEstimate, N,
#   tolerance = .25
# )

expect_true(
  (population_ztoipoisson$confidenceInterval[1, 1] < N) &
  (N < population_ztoipoisson$confidenceInterval[1, 2])
)

expect_true(
  (population_ztoipoisson$confidenceInterval[2, 1] < N) &
  (N < population_ztoipoisson$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_ztoipoisson,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztoipoisson,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztoipoisson,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_silent(
  predict(
    fit_ztoipoisson,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_ztoipoisson,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztoipoisson,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztoipoisson,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztoipoisson,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_ztoipoisson,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)
