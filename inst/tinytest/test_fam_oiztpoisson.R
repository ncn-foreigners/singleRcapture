source("helper_families.R")
if (!developer_tests_enabled()) exit_file("developer tests disabled")

inflated_test_data <- load_inflated_test_data()
N <- inflated_test_data$N
x1 <- inflated_test_data$x1
test_inflated <- inflated_test_data$test_inflated

# oizt ####

# poisson
coefficients_oiztpoisson <- c(.6, -.3, -.5, .1)
#eta <- cbind(coefficients_oiztpoisson[1] + coefficients_oiztpoisson[2] * x1, coefficients_oiztpoisson[3] + coefficients_oiztpoisson[4] * x1)

model_oiztpoisson <- oiztpoisson(
  lambdaLink = "log",
  omegaLink  = "logit",
)

#y <- model_oiztpoisson$simulate(n = N, eta = eta, lower = -1)
observed_oiztpoisson <- test_inflated$oiztpoisson1

data_oiztpoisson <- data.frame(
  x = x1[observed_oiztpoisson > 0],
  y = observed_oiztpoisson[observed_oiztpoisson > 0]
)

fit_oiztpoisson <- estimatePopsize(
  formula = y ~ x,
  model   = model_oiztpoisson,
  data    = data_oiztpoisson,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_oiztpoisson))

expect_silent(print(summary(fit_oiztpoisson)))

# expect_equivalent(
#   coef(fit_oiztpoisson),
#   coefficients_oiztpoisson,
#   tolerance = .25
# )

expect_silent(
  population_oiztpoisson <- popSizeEst(fit_oiztpoisson)
)

# expect_equivalent(
#   population_oiztpoisson$pointEstimate,
#   N,
#   tolerance = .2
# )

# expect_true(
#   (population_oiztpoisson$confidenceInterval[1, 1] < N) &
#   (N < population_oiztpoisson$confidenceInterval[1, 2])
# )

# expect_true(
#   (population_oiztpoisson$confidenceInterval[2, 1] < N) &
#   (N < population_oiztpoisson$confidenceInterval[2, 2])
# )

expect_true(
  all(predict(
    fit_oiztpoisson,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_oiztpoisson,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_oiztpoisson,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_oiztpoisson,
    type = "contr"
  )),
  fit_oiztpoisson$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_oiztpoisson,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_oiztpoisson,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_oiztpoisson,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_oiztpoisson,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_oiztpoisson,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_oiztpoisson,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)


# same with different link
model_oiztpoisson <- oiztpoisson(
  lambdaLink = "neglog",
  omegaLink  = "cloglog",
)

#y <- model_oiztpoisson$simulate(n = N, eta = eta, lower = -1)
observed_oiztpoisson <- test_inflated$oiztpoisson2

data_oiztpoisson <- data.frame(
  x = x1[observed_oiztpoisson > 0],
  y = observed_oiztpoisson[observed_oiztpoisson > 0]
)

fit_oiztpoisson <- estimatePopsize(
  formula = y ~ x,
  model   = model_oiztpoisson,
  data    = data_oiztpoisson,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_oiztpoisson))

expect_silent(print(summary(fit_oiztpoisson)))

# expect_equivalent(
#   coef(fit_oiztpoisson),
#   coefficients_oiztpoisson,
#   tolerance = .5
# )

expect_silent(
  population_oiztpoisson <- popSizeEst(fit_oiztpoisson)
)

# expect_equivalent(
#   population_oiztpoisson$pointEstimate,
#   N,
#   tolerance = .2
# )

expect_true(
  (population_oiztpoisson$confidenceInterval[1, 1] < N) &
  (N < population_oiztpoisson$confidenceInterval[1, 2])
)

# expect_true(
#   (population_oiztpoisson$confidenceInterval[2, 1] < N) &
#   (N < population_oiztpoisson$confidenceInterval[2, 2])
# )

expect_true(
  all(predict(
    fit_oiztpoisson,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_oiztpoisson,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_oiztpoisson,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_oiztpoisson,
    type = "contr"
  )),
  fit_oiztpoisson$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_oiztpoisson,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_oiztpoisson,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_oiztpoisson,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_oiztpoisson,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_oiztpoisson,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_oiztpoisson,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)
