source("helper_families.R")
if (!developer_tests_enabled()) exit_file("developer tests disabled")

inflated_test_data <- load_inflated_test_data()
N <- inflated_test_data$N
x1 <- inflated_test_data$x1
test_inflated <- inflated_test_data$test_inflated

# geometric
coefficients_ztHurdlegeom <- c(.6, .2, -3, .4)
#eta <- cbind(coefficients_ztHurdlegeom[1] + coefficients_ztHurdlegeom[2] * x1, coefficients_ztHurdlegeom[3] + coefficients_ztHurdlegeom[4] * x1)

model_ztHurdlegeom <- ztHurdlegeom(
  lambdaLink = "log",
  omegaLink  = "logit",
)

#y <- model_ztHurdlegeom$simulate(n = N, eta = eta, lower = -1)
observed_ztHurdlegeom <- test_inflated$ztHurdlegeom1

data_ztHurdlegeom <- data.frame(
  x = x1[observed_ztHurdlegeom > 0],
  y = observed_ztHurdlegeom[observed_ztHurdlegeom > 0]
)

fit_ztHurdlegeom <- estimatePopsize(
  formula = y ~ x,
  model   = model_ztHurdlegeom,
  data    = data_ztHurdlegeom,
  method = "IRLS",
  controlModel = controlModel(
    piFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_ztHurdlegeom))

expect_silent(print(summary(fit_ztHurdlegeom)))

expect_equivalent(
  coef(fit_ztHurdlegeom),
  coefficients_ztHurdlegeom,
  tolerance = .4
)

expect_silent(
  population_ztHurdlegeom <- popSizeEst(fit_ztHurdlegeom)
)

expect_equivalent(
  population_ztHurdlegeom$pointEstimate,
  N,
  tolerance = .2
)

expect_true(
  (population_ztHurdlegeom$confidenceInterval[1, 1] < N) &
  (N < population_ztHurdlegeom$confidenceInterval[1, 2])
)

expect_true(
  (population_ztHurdlegeom$confidenceInterval[2, 1] < N) &
  (N < population_ztHurdlegeom$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_ztHurdlegeom,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztHurdlegeom,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztHurdlegeom,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_ztHurdlegeom,
    type = "contr"
  )),
  fit_ztHurdlegeom$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztHurdlegeom,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_ztHurdlegeom,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztHurdlegeom,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztHurdlegeom,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztHurdlegeom,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_ztHurdlegeom,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)

# same for different link
coefficients_ztHurdlegeom <- c(.6, -.3, -1, .4)
#eta <- cbind(coefficients_ztHurdlegeom[1] + coefficients_ztHurdlegeom[2] * x1, coefficients_ztHurdlegeom[3] + coefficients_ztHurdlegeom[4] * x1)


model_ztHurdlegeom <- ztHurdlegeom(
  lambdaLink = "neglog",
  piLink = "probit",
)

#y <- model_ztHurdlegeom$simulate(n = N, eta = eta, lower = -1)
observed_ztHurdlegeom <- test_inflated$ztHurdlegeom2

data_ztHurdlegeom <- data.frame(
  x = x1[observed_ztHurdlegeom > 0],
  y = observed_ztHurdlegeom[observed_ztHurdlegeom > 0]
)

fit_ztHurdlegeom <- estimatePopsize(
  formula = y ~ x,
  model   = model_ztHurdlegeom,
  data    = data_ztHurdlegeom,
  method = "IRLS",
  controlModel = controlModel(
    piFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_ztHurdlegeom))

expect_silent(print(summary(fit_ztHurdlegeom)))

expect_equivalent(
  coef(fit_ztHurdlegeom),
  coefficients_ztHurdlegeom,
  tolerance = .65
)

expect_silent(
  population_ztHurdlegeom <- popSizeEst(fit_ztHurdlegeom)
)

expect_equivalent(
  population_ztHurdlegeom$pointEstimate,
  N,
  tolerance = .25
)

expect_true(
  (population_ztHurdlegeom$confidenceInterval[1, 1] < N) &
  (N < population_ztHurdlegeom$confidenceInterval[1, 2])
)

expect_true(
  (population_ztHurdlegeom$confidenceInterval[2, 1] < N) &
  (N < population_ztHurdlegeom$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_ztHurdlegeom,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztHurdlegeom,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztHurdlegeom,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_ztHurdlegeom,
    type = "contr"
  )),
  fit_ztHurdlegeom$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztHurdlegeom,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_ztHurdlegeom,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztHurdlegeom,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztHurdlegeom,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztHurdlegeom,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_ztHurdlegeom,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)
