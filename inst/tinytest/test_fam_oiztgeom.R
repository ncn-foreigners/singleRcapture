source("helper_families.R")
if (!developer_tests_enabled()) exit_file("developer tests disabled")

inflated_test_data <- load_inflated_test_data()
N <- inflated_test_data$N
x1 <- inflated_test_data$x1
test_inflated <- inflated_test_data$test_inflated

# geometric
coefficients_oiztgeom <- c(.6, .2, -1, .4)
#eta <- cbind(coefficients_oiztgeom[1] + coefficients_oiztgeom[2] * x1, coefficients_oiztgeom[3] + coefficients_oiztgeom[4] * x1)

model_oiztgeom <- oiztgeom(
  lambdaLink = "log",
  omegaLink  = "logit",
)

#y <- model_oiztgeom$simulate(n = N, eta = eta, lower = -1)
observed_oiztgeom <- test_inflated$oiztgeom1

data_oiztgeom <- data.frame(
  x = x1[observed_oiztgeom > 0],
  y = observed_oiztgeom[observed_oiztgeom > 0]
)

fit_oiztgeom <- estimatePopsize(
  formula = y ~ x,
  model   = model_oiztgeom,
  data    = data_oiztgeom,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_oiztgeom))

expect_silent(print(summary(fit_oiztgeom)))

expect_equivalent(
  coef(fit_oiztgeom),
  coefficients_oiztgeom,
  tolerance = .3
)

expect_silent(
  population_oiztgeom <- popSizeEst(fit_oiztgeom)
)

expect_equivalent(
  population_oiztgeom$pointEstimate,
  N,
  tolerance = .2
)


# expect_true(
#   (population_oiztgeom$confidenceInterval[1, 1] < N) &
#   (N < population_oiztgeom$confidenceInterval[1, 2])
# )

# expect_true(
#   (population_oiztgeom$confidenceInterval[2, 1] < N) &
#   (N < population_oiztgeom$confidenceInterval[2, 2])
# )

expect_true(
  all(predict(
    fit_oiztgeom,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_oiztgeom,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_oiztgeom,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_oiztgeom,
    type = "contr"
  )),
  fit_oiztgeom$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_oiztgeom,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_oiztgeom,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_oiztgeom,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_oiztgeom,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_oiztgeom,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_oiztgeom,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)

# same for different link
coefficients_oiztgeom <- c(.6, -.3, -1, .4)
#eta <- cbind(coefficients_oiztgeom[1] + coefficients_oiztgeom[2] * x1, coefficients_oiztgeom[3] + coefficients_oiztgeom[4] * x1)

model_oiztgeom <- oiztgeom(
  lambdaLink = "neglog",
  omegaLink  = "cloglog",
)

#y <- model_oiztgeom$simulate(n = N, eta = eta, lower = -1)
observed_oiztgeom <- test_inflated$oiztgeom2

data_oiztgeom <- data.frame(
  x = x1[observed_oiztgeom > 0],
  y = observed_oiztgeom[observed_oiztgeom > 0]
)

fit_oiztgeom <- estimatePopsize(
  formula = y ~ x,
  model   = model_oiztgeom,
  data    = data_oiztgeom,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_oiztgeom))

expect_silent(print(summary(fit_oiztgeom)))

expect_equivalent(
  coef(fit_oiztgeom),
  coefficients_oiztgeom,
  tolerance = .7
)

expect_silent(
  population_oiztgeom <- popSizeEst(fit_oiztgeom)
)

# expect_equivalent(
#   population_oiztgeom$pointEstimate,
#   N,
#   tolerance = .2
# )

# expect_true(
#   (population_oiztgeom$confidenceInterval[1, 1] < N) &
#   (N < population_oiztgeom$confidenceInterval[1, 2])
# )

# expect_true(
#   (population_oiztgeom$confidenceInterval[2, 1] < N) &
#   (N < population_oiztgeom$confidenceInterval[2, 2])
# )

expect_true(
  all(predict(
    fit_oiztgeom,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_oiztgeom,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_oiztgeom,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_oiztgeom,
    type = "contr"
  )),
  fit_oiztgeom$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_oiztgeom,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_oiztgeom,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_oiztgeom,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_oiztgeom,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_oiztgeom,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_oiztgeom,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)
