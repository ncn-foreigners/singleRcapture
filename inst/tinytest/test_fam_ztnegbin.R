source("helper_families.R")
zerotruncated_test_data <- load_zerotruncated_test_data()
eps <- zerotruncated_test_data$eps
N <- zerotruncated_test_data$N
test_zerotruncated <- zerotruncated_test_data$test_zerotruncated
x1 <- zerotruncated_test_data$x1
x2 <- zerotruncated_test_data$x2
x3 <- zerotruncated_test_data$x3

## ztnegbin ####
coefficients_ztnegbin <- c(.3, .2, -.3, .1,
          .2, .1, -.5, .3)
#eta <- cbind(coefficients_ztnegbin[1] + coefficients_ztnegbin[2] * x1 + coefficients_ztnegbin[3] * x2 + coefficients_ztnegbin[4] * x3,
#             coefficients_ztnegbin[5] + coefficients_ztnegbin[6] * x1 + coefficients_ztnegbin[7] * x2 + coefficients_ztnegbin[8] * x3)

model_ztnegbin <- ztnegbin(
  lambdaLink = "log",
  alphaLink = "log",
  eimStep = 10
)

#y <- model_ztnegbin$simulate(n = N, eta = eta, lower = -1)
observed_ztnegbin <- test_zerotruncated$ztnegbin1

data_ztnegbin <- data.frame(
  x3 = x3[observed_ztnegbin > 0],
  x2 = x2[observed_ztnegbin > 0],
  x1 = x1[observed_ztnegbin > 0],
  y = observed_ztnegbin[observed_ztnegbin > 0]
)

fit_ztnegbin <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = model_ztnegbin,
  data    = data_ztnegbin,
  method = "IRLS",
  controlModel = controlModel(
    alphaFormula = ~ x1 + x2 + x3
  )
)

expect_silent(
  population_ztnegbin <- popSizeEst(fit_ztnegbin)
)

expect_equivalent(
  population_ztnegbin$pointEstimate,
  N,
  tol = .3
)

expect_true(
  all(N < population_ztnegbin$confidenceInterval$upperBound) &
  all(N > population_ztnegbin$confidenceInterval$lowerBound)
)

expect_silent(
  strata_ztnegbin <- stratifyPopsize(fit_ztnegbin, ~ x1)
)

expect_equivalent(
  as.numeric(table(data_ztnegbin$x1))[as.numeric(substr(strata_ztnegbin$name, start = 5, stop = 5)) + 1],
  strata_ztnegbin$Observed
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(strata_ztnegbin$name, start = 5, stop = 5)) + 1],
  strata_ztnegbin$Estimated,
  tolerance = .5
)

# different link
coefficients_ztnegbin <- c(.1, .2, -.3, .1,
          .2, .1, -.1,   .3)

model_ztnegbin <- ztnegbin(
  lambdaLink = "log",
  alphaLink = "neglog",
  eimStep = 40
)

#y <- model_ztnegbin$simulate(n = N, eta = eta, lower = -1)
observed_ztnegbin <- test_zerotruncated$ztnegbin2

data_ztnegbin <- data.frame(
  x3 = x3[observed_ztnegbin > 0],
  x2 = x2[observed_ztnegbin > 0],
  x1 = x1[observed_ztnegbin > 0],
  y = observed_ztnegbin[observed_ztnegbin > 0]
)

fit_ztnegbin <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = model_ztnegbin,
  data    = data_ztnegbin,
  method = "IRLS",
  controlModel = controlModel(
    alphaFormula = ~ x1 + x2 + x3
  )
)

expect_true(
  max(abs(coef(fit_ztnegbin) - coefficients_ztnegbin)) < .9
)

expect_silent(
  population_ztnegbin <- popSizeEst(fit_ztnegbin)
)

expect_equivalent(
  N,
  population_ztnegbin$pointEstimate,
  tol = .1
)

expect_true(
  all(N < population_ztnegbin$confidenceInterval$upperBound) &
  all(N > population_ztnegbin$confidenceInterval$lowerBound)
)

expect_silent(
  strata_ztnegbin <- stratifyPopsize(fit_ztnegbin, ~ x1)
)

expect_equivalent(
  as.numeric(table(data_ztnegbin$x1))[as.numeric(substr(strata_ztnegbin$name, start = 5, stop = 5)) + 1],
  strata_ztnegbin$Observed,
  tolerance = eps
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(strata_ztnegbin$name, start = 5, stop = 5)) + 1],
  strata_ztnegbin$Estimated,
  tolerance = .1
)

marginal_ztnegbin <- summary(marginalFreq(fit_ztnegbin), dropl5 = "no")

expect_true(
  all(marginal_ztnegbin$Test$`P(>X^2)` > .05)
)

expect_true(
  all(predict(
    fit_ztnegbin,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztnegbin,
    type = "link",
    se.fit = TRUE
  )[,3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztnegbin,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    fit_ztnegbin,
    type = "contr"
  )) == fit_ztnegbin$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztnegbin,
    type = "popSize",
    se.fit = TRUE
  )
)

# different link
model_ztnegbin <- ztnegbin(
  lambdaLink = "neglog",
  alphaLink = "log",
  eimStep = 10
)

#y <- model_ztnegbin$simulate(n = N, eta = eta, lower = -1)
observed_ztnegbin <- test_zerotruncated$ztnegbin3

data_ztnegbin <- data.frame(
  x3 = x3[observed_ztnegbin > 0],
  x2 = x2[observed_ztnegbin > 0],
  x1 = x1[observed_ztnegbin > 0],
  y = observed_ztnegbin[observed_ztnegbin > 0]
)

fit_ztnegbin <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = model_ztnegbin,
  data    = data_ztnegbin,
  method = "IRLS",
  controlModel = controlModel(
    alphaFormula = ~ x1 + x2 + x3
  )
)

expect_silent(
  population_ztnegbin <- popSizeEst(fit_ztnegbin)
)

expect_equivalent(
  N,
  population_ztnegbin$pointEstimate,
  tol = .6
)

expect_true(
  (N < population_ztnegbin$confidenceInterval$upperBound[2]) &
  (N > population_ztnegbin$confidenceInterval$lowerBound[2])
)

expect_silent(
  strata_ztnegbin <- stratifyPopsize(fit_ztnegbin, ~ x1)
)

expect_equivalent(
  as.numeric(table(data_ztnegbin$x1))[as.numeric(substr(strata_ztnegbin$name, start = 5, stop = 5)) + 1],
  strata_ztnegbin$Observed,
  tolerance = eps
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(strata_ztnegbin$name, start = 5, stop = 5)) + 1],
  strata_ztnegbin$Estimated,
  tolerance = .6
)

expect_true(
  all(predict(
    fit_ztnegbin,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztnegbin,
    type = "link",
    se.fit = TRUE
  )[, 3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztnegbin,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)] >= 0)
)

expect_true(
  sum(predict(
    fit_ztnegbin,
    type = "contr"
  )) == fit_ztnegbin$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztnegbin,
    type = "popSize",
    se.fit = TRUE
  )
)

# different link
coefficients_ztnegbin <- c(.6, -.2, -1.25, -.1,
          .2, -.1, -.5,   .3)
#eta <- cbind(coefficients_ztnegbin[1] + coefficients_ztnegbin[2] * x1 + coefficients_ztnegbin[3] * x2 + coefficients_ztnegbin[4] * x3,
#             coefficients_ztnegbin[5] + coefficients_ztnegbin[6] * x1 + coefficients_ztnegbin[7] * x2 + coefficients_ztnegbin[8] * x3)

model_ztnegbin <- ztnegbin(
  lambdaLink = "log",
  alphaLink = "log",
  eimStep = 25
)

#y <- model_ztnegbin$simulate(n = N, eta = eta, lower = -1)
observed_ztnegbin <- test_zerotruncated$ztnegbin4

data_ztnegbin <- data.frame(
  x3 = x3[observed_ztnegbin > 0],
  x2 = x2[observed_ztnegbin > 0],
  x1 = x1[observed_ztnegbin > 0],
  y = observed_ztnegbin[observed_ztnegbin > 0]
)

fit_ztnegbin <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = model_ztnegbin,
  data    = data_ztnegbin,
  method = "IRLS",
  controlModel = controlModel(
    alphaFormula = ~ x1 + x2 + x3
  )
)

expect_equivalent(
  coefficients_ztnegbin,
  coef(fit_ztnegbin),
  tolerance = .4
)

expect_silent(
  population_ztnegbin <- popSizeEst(fit_ztnegbin)
)

expect_equivalent(
  N,
  population_ztnegbin$pointEstimate,
  tol = .15
)

expect_true(
  all(N < population_ztnegbin$confidenceInterval$upperBound) &
  all(N > population_ztnegbin$confidenceInterval$lowerBound)
)

expect_silent(
  strata_ztnegbin <- stratifyPopsize(fit_ztnegbin, ~ x1)
)

expect_equivalent(
  as.numeric(table(data_ztnegbin$x1))[as.numeric(substr(strata_ztnegbin$name, start = 5, stop = 5)) + 1],
  strata_ztnegbin$Observed,
  tolerance = eps
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(strata_ztnegbin$name, start = 5, stop = 5)) + 1],
  strata_ztnegbin$Estimated,
  tolerance = .15
)

marginal_ztnegbin <- summary(marginalFreq(fit_ztnegbin), dropl5 = "group")

expect_true(
  all(marginal_ztnegbin$Test$`P(>X^2)` > .05)
)

expect_true(
  all(predict(
    fit_ztnegbin,
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    fit_ztnegbin,
    type = "link",
    se.fit = TRUE
  )[, 3:4]> 0)
)

expect_true(
  all(predict(
    fit_ztnegbin,
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)] >= 0)
)

expect_true(
  sum(predict(
    fit_ztnegbin,
    type = "contr"
  )) == fit_ztnegbin$populationSize$pointEstimate
)

expect_silent(
  predict(
    fit_ztnegbin,
    type = "popSize",
    se.fit = TRUE
  )
)
