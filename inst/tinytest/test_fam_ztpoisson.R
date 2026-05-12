source("helper_families.R")
zerotruncated_test_data <- load_zerotruncated_test_data()
eps <- zerotruncated_test_data$eps
N <- zerotruncated_test_data$N
test_zerotruncated <- zerotruncated_test_data$test_zerotruncated
x1 <- zerotruncated_test_data$x1
x2 <- zerotruncated_test_data$x2
x3 <- zerotruncated_test_data$x3

## poisson ####
coefficients_ztpoisson <- c(.6, .2, -1.25, .1)
#eta <- cbind(coefficients_ztpoisson[1] + coefficients_ztpoisson[2] * x1 + coefficients_ztpoisson[3] * x2 + coefficients_ztpoisson[4] * x3)

model_ztpoisson <- ztpoisson(
  lambdaLink = "neglog"
)

#y <- model_ztpoisson$simulate(n = N, eta = eta, lower = -1)
observed_ztpoisson <- test_zerotruncated$ztpoisson1

data_ztpoisson <- data.frame(
  x3 = x3[observed_ztpoisson > 0],
  x2 = x2[observed_ztpoisson > 0],
  x1 = x1[observed_ztpoisson > 0],
  y = observed_ztpoisson[observed_ztpoisson > 0]
)

fit_ztpoisson <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = model_ztpoisson,
  data    = data_ztpoisson,
  method = "IRLS"
)

expect_equivalent(
  coef(fit_ztpoisson),
  coefficients_ztpoisson,
  tolerance = .25
)

expect_silent(
  population_ztpoisson <- popSizeEst(fit_ztpoisson)
)

expect_equivalent(
  N,
  population_ztpoisson$pointEstimate,
  tol = .15
)

expect_true(
  all(N < population_ztpoisson$confidenceInterval$upperBound) &
  all(N > population_ztpoisson$confidenceInterval$lowerBound)
)

expect_silent(
  strata_ztpoisson <- stratifyPopsize(fit_ztpoisson, ~ x1)
)

expect_equivalent(
  as.numeric(table(data_ztpoisson$x1))[as.numeric(substr(strata_ztpoisson$name, start = 5, stop = 5)) + 1],
  strata_ztpoisson$Observed,
  tolerance = eps
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(strata_ztpoisson$name, start = 5, stop = 5)) + 1],
  strata_ztpoisson$Estimated,
  tolerance = .1
)

marginal_ztpoisson <- summary(marginalFreq(fit_ztpoisson), dropl5 = "group", df = 1)

expect_true(
  all(marginal_ztpoisson$Test$`P(>X^2)` > .01)
)

expect_true(
  all(predict(
    fit_ztpoisson,
    type = "response",
    newdata = test_zerotruncated[, paste0("x", 1:3)],
    se.fit = TRUE
  )[, 2] > 0)
)

expect_true(
  all(predict(
    fit_ztpoisson,
    type = "link",
    newdata = test_zerotruncated[, paste0("x", 1:3)],
    se.fit = TRUE
  )[,2]> 0)
)

expect_true(
  all(predict(
    fit_ztpoisson,
    type = "mean",
    newdata = test_zerotruncated[, paste0("x", 1:3)],
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztpoisson,
    type = "contr",
    newdata = test_zerotruncated[, paste0("x", 1:3)]
  )) > N
)

expect_silent(
  predict(
    fit_ztpoisson,
    type = "popSize",
    newdata = test_zerotruncated[, paste0("x", 1:3)],
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    fit_ztpoisson,
    type = "response",
    se.fit = TRUE
  )[, 2] > 0)
)

expect_true(
  all(predict(
    fit_ztpoisson,
    type = "link",
    se.fit = TRUE
  )[,2]> 0)
)

expect_true(
  all(predict(
    fit_ztpoisson,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztpoisson,
    type = "contr"
  )) == fit_ztpoisson$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztpoisson,
    type = "popSize",
    se.fit = TRUE
  )
)

# different link

model_ztpoisson <- ztpoisson(
  lambdaLink = "log"
)

#y <- model_ztpoisson$simulate(n = N, eta = eta, lower = -1)
observed_ztpoisson <- test_zerotruncated$ztpoisson2

data_ztpoisson <- data.frame(
  x3 = x3[observed_ztpoisson > 0],
  x2 = x2[observed_ztpoisson > 0],
  x1 = x1[observed_ztpoisson > 0],
  y = observed_ztpoisson[observed_ztpoisson > 0]
)

fit_ztpoisson <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = model_ztpoisson,
  data    = data_ztpoisson,
  method = "IRLS"
)

expect_equivalent(
  coef(fit_ztpoisson),
  coefficients_ztpoisson,
  tolerance = .15
)

expect_silent(
  population_ztpoisson <- popSizeEst(fit_ztpoisson)
)

expect_equivalent(
  N,
  population_ztpoisson$pointEstimate,
  tol = .1
)

expect_true(
  all(N < population_ztpoisson$confidenceInterval$upperBound) &
  all(N > population_ztpoisson$confidenceInterval$lowerBound)
)

expect_silent(
  strata_ztpoisson <- stratifyPopsize(fit_ztpoisson, ~ x1)
)

expect_equivalent(
  as.numeric(table(data_ztpoisson$x1))[as.numeric(substr(strata_ztpoisson$name, start = 5, stop = 5)) + 1],
  strata_ztpoisson$Observed,
  tolerance = eps
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(strata_ztpoisson$name, start = 5, stop = 5)) + 1],
  strata_ztpoisson$Estimated,
  tolerance = .1
)

marginal_ztpoisson <- summary(marginalFreq(fit_ztpoisson), dropl5 = "group")

expect_true(
  all(marginal_ztpoisson$Test$`P(>X^2)` > .01)
)

expect_true(
  all(predict(
    fit_ztpoisson,
    type = "response",
    se.fit = TRUE
  )[, 2] > 0)
)

expect_true(
  all(predict(
    fit_ztpoisson,
    type = "link",
    se.fit = TRUE
  )[,2]> 0)
)

expect_true(
  all(predict(
    fit_ztpoisson,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztpoisson,
    type = "contr"
  )) == fit_ztpoisson$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztpoisson,
    type = "popSize",
    se.fit = TRUE
  )
)

## bootstrap

fit_ztpoisson <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = model_ztpoisson,
  data    = data_ztpoisson,
  method = "IRLS",
  popVar = "bootstrap",
  controlPopVar = controlPopVar(
    B = 250,
    bootType = "parametric"
  )
)

expect_silent(
  plot(fit_ztpoisson, "bootHist")
)

expect_silent(
  print(summary(fit_ztpoisson, correlation = TRUE))
)

expect_silent(
  population_ztpoisson <- popSizeEst(fit_ztpoisson)
)

expect_equivalent(
  N,
  population_ztpoisson$pointEstimate,
  tol = .1
)

expect_true(
  all(N < population_ztpoisson$confidenceInterval[2]) &
  all(N > population_ztpoisson$confidenceInterval[1])
)

expect_true(
  shapiro.test(population_ztpoisson$boot)$p.value > .05
)

fit_ztpoisson <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = model_ztpoisson,
  data    = data_ztpoisson,
  method = "IRLS",
  popVar = "bootstrap",
  controlPopVar = controlPopVar(
    B = 250,
    bootType = "semiparametric"
  )
)

expect_silent(
  plot(fit_ztpoisson, "bootHist")
)

expect_silent(
  print(summary(fit_ztpoisson, correlation = TRUE))
)

expect_silent(
  population_ztpoisson <- popSizeEst(fit_ztpoisson)
)

expect_equivalent(
  N,
  population_ztpoisson$pointEstimate,
  tol = .1
)

expect_true(
  all(N < population_ztpoisson$confidenceInterval[2]) &
  all(N > population_ztpoisson$confidenceInterval[1])
)

expect_true(
  shapiro.test(population_ztpoisson$boot)$p.value > .05
)

expect_true(
  all(predict(
    fit_ztpoisson,
    type = "response",
    se.fit = TRUE
  )[, 2] > 0)
)

expect_true(
  all(predict(
    fit_ztpoisson,
    type = "link",
    se.fit = TRUE
  )[,2]> 0)
)

expect_true(
  all(predict(
    fit_ztpoisson,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztpoisson,
    type = "contr"
  )) == fit_ztpoisson$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztpoisson,
    type = "popSize",
    se.fit = TRUE
  )
)

fit_ztpoisson <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = model_ztpoisson,
  data    = data_ztpoisson,
  method = "IRLS",
  popVar = "bootstrap",
  controlPopVar = controlPopVar(
    B = 250,
    bootType = "nonparametric"
  )
)

expect_silent(
  plot(fit_ztpoisson, "bootHist")
)

expect_silent(
  print(summary(fit_ztpoisson, correlation = TRUE))
)

expect_silent(
  population_ztpoisson <- popSizeEst(fit_ztpoisson)
)

expect_equivalent(
  N,
  population_ztpoisson$pointEstimate,
  tol = .1
)
