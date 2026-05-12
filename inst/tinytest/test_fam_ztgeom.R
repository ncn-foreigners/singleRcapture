source("helper_families.R")
zerotruncated_test_data <- load_zerotruncated_test_data()
eps <- zerotruncated_test_data$eps
N <- zerotruncated_test_data$N
test_zerotruncated <- zerotruncated_test_data$test_zerotruncated
x1 <- zerotruncated_test_data$x1
x2 <- zerotruncated_test_data$x2
x3 <- zerotruncated_test_data$x3

## geometric ####
coefficients_ztgeom <- c(.6, .2, -1.25, .1)
#eta <- cbind(coefficients_ztgeom[1] + coefficients_ztgeom[2] * x1 + coefficients_ztgeom[3] * x2 + coefficients_ztgeom[4] * x3)

model_ztgeom <- ztgeom(
  lambdaLink = "neglog"
)

#y <- model_ztgeom$simulate(n = N, eta = eta, lower = -1)
observed_ztgeom <- test_zerotruncated$ztgeom1

data_ztgeom <- data.frame(
  x3 = x3[observed_ztgeom > 0],
  x2 = x2[observed_ztgeom > 0],
  x1 = x1[observed_ztgeom > 0],
  y = observed_ztgeom[observed_ztgeom > 0]
)

fit_ztgeom <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = model_ztgeom,
  data    = data_ztgeom,
  method = "IRLS"
)

expect_equivalent(
  coef(fit_ztgeom),
  coefficients_ztgeom,
  tolerance = .25
)

expect_silent(
  population_ztgeom <- popSizeEst(fit_ztgeom)
)

expect_equivalent(
  N,
  population_ztgeom$pointEstimate,
  tol = .1
)

expect_true(
  all(N < population_ztgeom$confidenceInterval$upperBound) &
  all(N > population_ztgeom$confidenceInterval$lowerBound)
)

expect_silent(
  strata_ztgeom <- stratifyPopsize(fit_ztgeom, ~ x1)
)

expect_equivalent(
  as.numeric(table(data_ztgeom$x1))[as.numeric(substr(strata_ztgeom$name, start = 5, stop = 5)) + 1],
  strata_ztgeom$Observed,
  tolerance = eps
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(strata_ztgeom$name, start = 5, stop = 5)) + 1],
  strata_ztgeom$Estimated,
  tolerance = .2
)

marginal_ztgeom <- summary(marginalFreq(fit_ztgeom), dropl5 = "group")

expect_true(
  all(marginal_ztgeom$Test$`P(>X^2)` > .03)
)

expect_true(
  all(predict(
    fit_ztgeom,
    type = "response",
    se.fit = TRUE
  )[, 2] > 0)
)

expect_true(
  all(predict(
    fit_ztgeom,
    type = "link",
    se.fit = TRUE
  )[,2]> 0)
)

expect_true(
  all(predict(
    fit_ztgeom,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztgeom,
    type = "contr"
  )) == fit_ztgeom$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztgeom,
    type = "popSize",
    se.fit = TRUE
  )
)

# different link

model_ztgeom <- ztgeom(
  lambdaLink = "log"
)

#y <- model_ztgeom$simulate(n = N, eta = eta, lower = -1)
observed_ztgeom <- test_zerotruncated$ztgeom2

data_ztgeom <- data.frame(
  x3 = x3[observed_ztgeom > 0],
  x2 = x2[observed_ztgeom > 0],
  x1 = x1[observed_ztgeom > 0],
  y = observed_ztgeom[observed_ztgeom > 0]
)

fit_ztgeom <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = model_ztgeom,
  data    = data_ztgeom,
  method = "IRLS"
)

expect_equivalent(
  coef(fit_ztgeom),
  coefficients_ztgeom,
  tolerance = .2
)

expect_silent(
  population_ztgeom <- popSizeEst(fit_ztgeom)
)

expect_equivalent(
  N,
  population_ztgeom$pointEstimate,
  tol = .05
)

expect_true(
  all(N < population_ztgeom$confidenceInterval$upperBound) &
  all(N > population_ztgeom$confidenceInterval$lowerBound)
)

expect_silent(
  strata_ztgeom <- stratifyPopsize(fit_ztgeom, ~ x1)
)

expect_equivalent(
  as.numeric(table(data_ztgeom$x1))[as.numeric(substr(strata_ztgeom$name, start = 5, stop = 5)) + 1],
  strata_ztgeom$Observed,
  tolerance = .1
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(strata_ztgeom$name, start = 5, stop = 5)) + 1],
  strata_ztgeom$Estimated,
  tolerance = .1
)

marginal_ztgeom <- summary(marginalFreq(fit_ztgeom), dropl5 = "group")

expect_true(
  all(marginal_ztgeom$Test$`P(>X^2)` > .05)
)

expect_true(
  all(predict(
    fit_ztgeom,
    type = "response",
    se.fit = TRUE
  )[, 2] > 0)
)

expect_true(
  all(predict(
    fit_ztgeom,
    type = "link",
    se.fit = TRUE
  )[,2]> 0)
)

expect_true(
  all(predict(
    fit_ztgeom,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztgeom,
    type = "contr"
  )) == fit_ztgeom$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztgeom,
    type = "popSize",
    se.fit = TRUE
  )
)
