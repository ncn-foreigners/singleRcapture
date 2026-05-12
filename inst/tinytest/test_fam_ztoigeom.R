source("helper_families.R")
if (!developer_tests_enabled()) exit_file("developer tests disabled")

inflated_test_data <- load_inflated_test_data()
N <- inflated_test_data$N
x1 <- inflated_test_data$x1
test_inflated <- inflated_test_data$test_inflated

# ztoigeom
coefficients_ztoigeom <- c(.6, .2, -1.25, .1)
#eta <- cbind(coefficients_ztoigeom[1] + coefficients_ztoigeom[2] * x1, coefficients_ztoigeom[3] + coefficients_ztoigeom[4] * x1)

model_ztoigeom <- ztoigeom(
  lambdaLink = "log",
  omegaLink  = "logit",
)

#y <- model_ztoigeom$simulate(n = N, eta = eta, lower = -1)
observed_ztoigeom <- test_inflated$ztoigeom1

data_ztoigeom <- data.frame(
  x = x1[observed_ztoigeom > 0],
  y = observed_ztoigeom[observed_ztoigeom > 0]
)

fit_ztoigeom <- estimatePopsize(
  formula = y ~ x,
  model   = model_ztoigeom,
  data    = data_ztoigeom,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_ztoigeom))

expect_silent(print(summary(fit_ztoigeom)))

expect_equivalent(
  coef(fit_ztoigeom),
  coefficients_ztoigeom,
  tolerance = .7
)

expect_silent(
  population_ztoigeom <- popSizeEst(fit_ztoigeom)
)

expect_equivalent(
  population_ztoigeom$pointEstimate,
  N,
  tolerance = .2
)

expect_true(
  (population_ztoigeom$confidenceInterval[1, 1] < N) &
  (N < population_ztoigeom$confidenceInterval[1, 2])
)

expect_true(
  (population_ztoigeom$confidenceInterval[2, 1] < N) &
  (N < population_ztoigeom$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_ztoigeom,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztoigeom,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztoigeom,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_ztoigeom,
    type = "contr"
  )),
  fit_ztoigeom$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztoigeom,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_ztoigeom,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztoigeom,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztoigeom,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztoigeom,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_ztoigeom,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)

# same for different link
model_ztoigeom <- ztoigeom(
  lambdaLink = "neglog",
  omegaLink  = "cloglog",
)

#y <- model_ztoigeom$simulate(n = N, eta = eta, lower = -1)
observed_ztoigeom <- test_inflated$ztoigeom2

data_ztoigeom <- data.frame(
  x = x1[observed_ztoigeom > 0],
  y = observed_ztoigeom[observed_ztoigeom > 0]
)

fit_ztoigeom <- estimatePopsize(
  formula = y ~ x,
  model   = model_ztoigeom,
  data    = data_ztoigeom,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_ztoigeom))

expect_silent(print(summary(fit_ztoigeom)))

expect_silent(
  population_ztoigeom <- popSizeEst(fit_ztoigeom)
)

expect_equivalent(
  population_ztoigeom$pointEstimate,
  N,
  tolerance = .2
)

expect_true(
  (population_ztoigeom$confidenceInterval[1, 1] < N) &
  (N < population_ztoigeom$confidenceInterval[1, 2])
)

expect_true(
  (population_ztoigeom$confidenceInterval[2, 1] < N) &
  (N < population_ztoigeom$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_ztoigeom,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztoigeom,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztoigeom,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_ztoigeom,
    type = "contr"
  )),
  fit_ztoigeom$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztoigeom,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_ztoigeom,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztoigeom,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztoigeom,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztoigeom,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_ztoigeom,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)
