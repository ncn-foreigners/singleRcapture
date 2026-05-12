source("helper_families.R")
if (!developer_tests_enabled()) exit_file("developer tests disabled")

inflated_test_data <- load_inflated_test_data()
N <- inflated_test_data$N
x1 <- inflated_test_data$x1
test_inflated <- inflated_test_data$test_inflated

# geometric
coefficients_Hurdleztgeom <- c(.6, .2, -3, .4)
#eta <- cbind(coefficients_Hurdleztgeom[1] + coefficients_Hurdleztgeom[2] * x1, coefficients_Hurdleztgeom[3] + coefficients_Hurdleztgeom[4] * x1)

model_Hurdleztgeom <- Hurdleztgeom(
  lambdaLink = "log",
  omegaLink  = "logit",
)

#y <- model_Hurdleztgeom$simulate(n = N, eta = eta, lower = -1)
observed_Hurdleztgeom <- test_inflated$Hurdleztgeom1

data_Hurdleztgeom <- data.frame(
  x = x1[observed_Hurdleztgeom > 0],
  y = observed_Hurdleztgeom[observed_Hurdleztgeom > 0]
)

fit_Hurdleztgeom <- estimatePopsize(
  formula = y ~ x,
  model   = model_Hurdleztgeom,
  data    = data_Hurdleztgeom,
  method = "IRLS",
  controlModel = controlModel(
    piFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_Hurdleztgeom))

expect_silent(print(summary(fit_Hurdleztgeom)))

expect_equivalent(
  coef(fit_Hurdleztgeom),
  coefficients_Hurdleztgeom,
  tolerance = .2
)

expect_silent(
  population_Hurdleztgeom <- popSizeEst(fit_Hurdleztgeom)
)

expect_equivalent(
  population_Hurdleztgeom$pointEstimate,
  N,
  tolerance = .2
)

expect_true(
  (population_Hurdleztgeom$confidenceInterval[1, 1] < N) &
  (N < population_Hurdleztgeom$confidenceInterval[1, 2])
)

expect_true(
  (population_Hurdleztgeom$confidenceInterval[2, 1] < N) &
  (N < population_Hurdleztgeom$confidenceInterval[2, 2])
)

# same for different link
coefficients_Hurdleztgeom <- c(.3, -.7, -.25, .3)
#eta <- cbind(coefficients_Hurdleztgeom[1] + coefficients_Hurdleztgeom[2] * x1, coefficients_Hurdleztgeom[3] + coefficients_Hurdleztgeom[4] * x1)

model_Hurdleztgeom <- Hurdleztgeom(
  lambdaLink = "neglog",
  piLink = "probit",
)

#y <- model_Hurdleztgeom$simulate(n = N, eta = eta, lower = -1)
observed_Hurdleztgeom <- test_inflated$Hurdleztgeom2

data_Hurdleztgeom <- data.frame(
  x = x1[observed_Hurdleztgeom > 0],
  y = observed_Hurdleztgeom[observed_Hurdleztgeom > 0]
)

fit_Hurdleztgeom <- estimatePopsize(
  formula = y ~ x,
  model   = model_Hurdleztgeom,
  data    = data_Hurdleztgeom,
  method = "IRLS",
  controlModel = controlModel(
    piFormula = ~ x
  ),
  controlMethod = controlMethod(silent = TRUE)
)

expect_silent(summary(fit_Hurdleztgeom))

expect_silent(print(summary(fit_Hurdleztgeom)))

expect_equivalent(
  coef(fit_Hurdleztgeom),
  coefficients_Hurdleztgeom,
  tolerance = .4
)

expect_silent(
  population_Hurdleztgeom <- popSizeEst(fit_Hurdleztgeom)
)

expect_equivalent(
  population_Hurdleztgeom$pointEstimate,
  N,
  tolerance = .2
)


expect_true(
  (population_Hurdleztgeom$confidenceInterval[1, 1] < N) &
  (N < population_Hurdleztgeom$confidenceInterval[1, 2])
)

expect_true(
  (population_Hurdleztgeom$confidenceInterval[2, 1] < N) &
  (N < population_Hurdleztgeom$confidenceInterval[2, 2])
)

expect_true(
  all(predict(
    fit_Hurdleztgeom,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_Hurdleztgeom,
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_Hurdleztgeom,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_equivalent(
  sum(predict(
    fit_Hurdleztgeom,
    type = "contr"
  )),
  fit_Hurdleztgeom$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_Hurdleztgeom,
    type = "popSize",
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_Hurdleztgeom,
    newdata = data.frame(x = x1),
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_Hurdleztgeom,
    newdata = data.frame(x = x1),
    type = "link",
    se.fit = TRUE
  )[,3:4]> 0)
)

expect_true(
  all(predict(
    fit_Hurdleztgeom,
    newdata = data.frame(x = x1),
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_Hurdleztgeom,
    newdata = data.frame(x = x1),
    type = "contr"
  )) > N
)

expect_silent(
  predict(
    fit_Hurdleztgeom,
    newdata = data.frame(x = x1),
    type = "popSize",
    se.fit = TRUE
  )
)
