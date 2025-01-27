# testing for zero truncated models

eps <- if (capabilities("long.double")) sqrt(.Machine$double.eps) else 0.01

set.seed(123)

N <- 1000
#x1 <- rbinom(n = N, size = 3, prob = .4)
#x2 <- runif(n = N)
#x3 <- rnorm(n = N, mean = 3, sd = 4)

testTruncated <- read.csv("test_zerotruncated.csv")
#testTruncated <- read.csv("inst/tinytest/test_zerotruncated.csv")
x1 <- testTruncated$x1
x2 <- testTruncated$x2
x3 <- testTruncated$x3

## poisson ####
beta <- c(.6, .2, -1.25, .1)
#eta <- cbind(beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x3)

fn <- ztpoisson(
  lambdaLink = "neglog"
)

#y <- fn$simulate(n = N, eta = eta, lower = -1)
y <- testTruncated$ztpoisson1

df <- data.frame(
  x3 = x3[y > 0],
  x2 = x2[y > 0],
  x1 = x1[y > 0],
  y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = fn,
  data    = df,
  method = "IRLS"
)

expect_equivalent(
  coef(M1),
  beta,
  tolerance = .25
)

expect_silent(
  pop <- popSizeEst(M1)
)

expect_equivalent(
  N,
  pop$pointEstimate,
  tol = .15
)

expect_true(
  all(N < pop$confidenceInterval$upperBound) &
  all(N > pop$confidenceInterval$lowerBound)
)

expect_silent(
  stra <- stratifyPopsize(M1, ~ x1)
)

expect_equivalent(
  as.numeric(table(df$x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Observed,
  tolerance = eps
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Estimated,
  tolerance = .1
)

AA <- summary(marginalFreq(M1), dropl5 = "group", df = 1)

expect_true(
  all(AA$Test$`P(>X^2)` > .01)
)

expect_true(
  all(predict(
    M1, 
    type = "response",
    newdata = testTruncated[, paste0("x", 1:3)],
    se.fit = TRUE
  )[, 2] > 0)
)

expect_true(
  all(predict(
    M1, 
    type = "link",
    newdata = testTruncated[, paste0("x", 1:3)],
    se.fit = TRUE
  )[,2]> 0)
)

expect_true(
  all(predict(
    M1, 
    type = "mean",
    newdata = testTruncated[, paste0("x", 1:3)],
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    M1, 
    type = "contr",
    newdata = testTruncated[, paste0("x", 1:3)]
  )) > N
)

expect_silent(
  predict(
    M1, 
    type = "popSize",
    newdata = testTruncated[, paste0("x", 1:3)],
    se.fit = TRUE
  )
)

expect_true(
  all(predict(
    M1, 
    type = "response",
    se.fit = TRUE
  )[, 2] > 0)
)

expect_true(
  all(predict(
    M1, 
    type = "link",
    se.fit = TRUE
  )[,2]> 0)
)

expect_true(
  all(predict(
    M1, 
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    M1, 
    type = "contr"
  )) == M1$populationSize$pointEstimate
)

expect_silent(
  predict(
    M1, 
    type = "popSize",
    se.fit = TRUE
  )
)

# different link

fn <- ztpoisson(
  lambdaLink = "log"
)

#y <- fn$simulate(n = N, eta = eta, lower = -1)
y <- testTruncated$ztpoisson2

df <- data.frame(
  x3 = x3[y > 0],
  x2 = x2[y > 0],
  x1 = x1[y > 0],
  y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = fn,
  data    = df,
  method = "IRLS"
)

expect_equivalent(
  coef(M1),
  beta,
  tolerance = .15
)

expect_silent(
  pop <- popSizeEst(M1)
)

expect_equivalent(
  N,
  pop$pointEstimate,
  tol = .1
)

expect_true(
  all(N < pop$confidenceInterval$upperBound) &
  all(N > pop$confidenceInterval$lowerBound)
)

expect_silent(
  stra <- stratifyPopsize(M1, ~ x1)
)

expect_equivalent(
  as.numeric(table(df$x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Observed,
  tolerance = eps
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Estimated,
  tolerance = .1
)

AA <- summary(marginalFreq(M1), dropl5 = "group")

expect_true(
  all(AA$Test$`P(>X^2)` > .01)
)

expect_true(
  all(predict(
    M1, 
    type = "response",
    se.fit = TRUE
  )[, 2] > 0)
)

expect_true(
  all(predict(
    M1, 
    type = "link",
    se.fit = TRUE
  )[,2]> 0)
)

expect_true(
  all(predict(
    M1, 
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    M1, 
    type = "contr"
  )) == M1$populationSize$pointEstimate
)

expect_silent(
  predict(
    M1, 
    type = "popSize",
    se.fit = TRUE
  )
)

## bootstrap

M1 <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = fn,
  data    = df,
  method = "IRLS",
  popVar = "bootstrap",
  controlPopVar = controlPopVar(
    B = 250,
    bootType = "parametric"
  )
)

expect_silent(
  plot(M1, "bootHist")
)

expect_silent(
  print(summary(M1, correlation = TRUE))
)

expect_silent(
  pop <- popSizeEst(M1)
)

expect_equivalent(
  N,
  pop$pointEstimate,
  tol = .1
)

expect_true(
  all(N < pop$confidenceInterval[2]) &
  all(N > pop$confidenceInterval[1])
)

expect_true(
  shapiro.test(pop$boot)$p.value > .05
)

M1 <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = fn,
  data    = df,
  method = "IRLS",
  popVar = "bootstrap",
  controlPopVar = controlPopVar(
    B = 250,
    bootType = "semiparametric"
  )
)

expect_silent(
  plot(M1, "bootHist")
)

expect_silent(
  print(summary(M1, correlation = TRUE))
)

expect_silent(
  pop <- popSizeEst(M1)
)

expect_equivalent(
  N,
  pop$pointEstimate,
  tol = .1
)

expect_true(
  all(N < pop$confidenceInterval[2]) &
  all(N > pop$confidenceInterval[1])
)

expect_true(
  shapiro.test(pop$boot)$p.value > .05
)

expect_true(
  all(predict(
    M1, 
    type = "response",
    se.fit = TRUE
  )[, 2] > 0)
)

expect_true(
  all(predict(
    M1, 
    type = "link",
    se.fit = TRUE
  )[,2]> 0)
)

expect_true(
  all(predict(
    M1, 
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    M1, 
    type = "contr"
  )) == M1$populationSize$pointEstimate
)

expect_silent(
  predict(
    M1, 
    type = "popSize",
    se.fit = TRUE
  )
)

M1 <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = fn,
  data    = df,
  method = "IRLS",
  popVar = "bootstrap",
  controlPopVar = controlPopVar(
    B = 250,
    bootType = "nonparametric"
  )
)

expect_silent(
  plot(M1, "bootHist")
)

expect_silent(
  print(summary(M1, correlation = TRUE))
)

expect_silent(
  pop <- popSizeEst(M1)
)

expect_equivalent(
  N,
  pop$pointEstimate,
  tol = .1
)

## geometric ####
beta <- c(.6, .2, -1.25, .1)
#eta <- cbind(beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x3)

fn <- ztgeom(
  lambdaLink = "neglog"
)

#y <- fn$simulate(n = N, eta = eta, lower = -1)
y <- testTruncated$ztgeom1

df <- data.frame(
  x3 = x3[y > 0],
  x2 = x2[y > 0],
  x1 = x1[y > 0],
  y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = fn,
  data    = df,
  method = "IRLS"
)

expect_equivalent(
  coef(M1),
  beta,
  tolerance = .25
)

expect_silent(
  pop <- popSizeEst(M1)
)

expect_equivalent(
  N,
  pop$pointEstimate,
  tol = .1
)

expect_true(
  all(N < pop$confidenceInterval$upperBound) &
  all(N > pop$confidenceInterval$lowerBound)
)

expect_silent(
  stra <- stratifyPopsize(M1, ~ x1)
)

expect_equivalent(
  as.numeric(table(df$x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Observed,
  tolerance = eps
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Estimated,
  tolerance = .2
)

AA <- summary(marginalFreq(M1), dropl5 = "group")

expect_true(
  all(AA$Test$`P(>X^2)` > .03)
)

expect_true(
  all(predict(
    M1, 
    type = "response",
    se.fit = TRUE
  )[, 2] > 0)
)

expect_true(
  all(predict(
    M1, 
    type = "link",
    se.fit = TRUE
  )[,2]> 0)
)

expect_true(
  all(predict(
    M1, 
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    M1, 
    type = "contr"
  )) == M1$populationSize$pointEstimate
)

expect_silent(
  predict(
    M1, 
    type = "popSize",
    se.fit = TRUE
  )
)

# different link

fn <- ztgeom(
  lambdaLink = "log"
)

#y <- fn$simulate(n = N, eta = eta, lower = -1)
y <- testTruncated$ztgeom2

df <- data.frame(
  x3 = x3[y > 0],
  x2 = x2[y > 0],
  x1 = x1[y > 0],
  y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = fn,
  data    = df,
  method = "IRLS"
)

expect_equivalent(
  coef(M1),
  beta,
  tolerance = .2
)

expect_silent(
  pop <- popSizeEst(M1)
)

expect_equivalent(
  N,
  pop$pointEstimate,
  tol = .05
)

expect_true(
  all(N < pop$confidenceInterval$upperBound) &
  all(N > pop$confidenceInterval$lowerBound)
)

expect_silent(
  stra <- stratifyPopsize(M1, ~ x1)
)

expect_equivalent(
  as.numeric(table(df$x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Observed,
  tolerance = .1
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Estimated,
  tolerance = .1
)

AA <- summary(marginalFreq(M1), dropl5 = "group")

expect_true(
  all(AA$Test$`P(>X^2)` > .05)
)

expect_true(
  all(predict(
    M1, 
    type = "response",
    se.fit = TRUE
  )[, 2] > 0)
)

expect_true(
  all(predict(
    M1, 
    type = "link",
    se.fit = TRUE
  )[,2]> 0)
)

expect_true(
  all(predict(
    M1, 
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    M1, 
    type = "contr"
  )) == M1$populationSize$pointEstimate
)

expect_silent(
  predict(
    M1, 
    type = "popSize",
    se.fit = TRUE
  )
)

## ztnegbin ####
beta <- c(.3, .2, -.3, .1,
          .2, .1, -.5, .3)
#eta <- cbind(beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x3,
#             beta[5] + beta[6] * x1 + beta[7] * x2 + beta[8] * x3)

fn <- ztnegbin(
  lambdaLink = "log",
  alphaLink = "log",
  eimStep = 10
)

#y <- fn$simulate(n = N, eta = eta, lower = -1)
y <- testTruncated$ztnegbin1

df <- data.frame(
  x3 = x3[y > 0],
  x2 = x2[y > 0],
  x1 = x1[y > 0],
  y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = fn,
  data    = df,
  method = "IRLS",
  controlModel = controlModel(
    alphaFormula = ~ x1 + x2 + x3
  )
)

expect_silent(
  pop <- popSizeEst(M1)
)

expect_equivalent(
  pop$pointEstimate,
  N,
  tol = .3
)

expect_true(
  all(N < pop$confidenceInterval$upperBound) &
  all(N > pop$confidenceInterval$lowerBound)
)

expect_silent(
  stra <- stratifyPopsize(M1, ~ x1)
)

expect_equivalent(
  as.numeric(table(df$x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Observed
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Estimated,
  tolerance = .5
)

# different link
beta <- c(.1, .2, -.3, .1,
          .2, .1, -.1,   .3)

fn <- ztnegbin(
  lambdaLink = "log",
  alphaLink = "neglog",
  eimStep = 40
)

#y <- fn$simulate(n = N, eta = eta, lower = -1)
y <- testTruncated$ztnegbin2

df <- data.frame(
  x3 = x3[y > 0],
  x2 = x2[y > 0],
  x1 = x1[y > 0],
  y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = fn,
  data    = df,
  method = "IRLS",
  controlModel = controlModel(
    alphaFormula = ~ x1 + x2 + x3
  )
)

expect_true(
  max(abs(coef(M1) - beta)) < .9
)

expect_silent(
  pop <- popSizeEst(M1)
)

expect_equivalent(
  N,
  pop$pointEstimate,
  tol = .1
)

expect_true(
  all(N < pop$confidenceInterval$upperBound) &
  all(N > pop$confidenceInterval$lowerBound)
)

expect_silent(
  stra <- stratifyPopsize(M1, ~ x1)
)

expect_equivalent(
  as.numeric(table(df$x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Observed,
  tolerance = eps
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Estimated,
  tolerance = .1
)

AA <- summary(marginalFreq(M1), dropl5 = "no")

expect_true(
  all(AA$Test$`P(>X^2)` > .05)
)

expect_true(
  all(predict(
    M1, 
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    M1, 
    type = "link",
    se.fit = TRUE
  )[,3:4] > 0)
)

expect_true(
  all(predict(
    M1, 
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)]> 0)
)

expect_true(
  sum(predict(
    M1, 
    type = "contr"
  )) == M1$populationSize$pointEstimate
)

expect_silent(
  predict(
    M1, 
    type = "popSize",
    se.fit = TRUE
  )
)

# different link
fn <- ztnegbin(
  lambdaLink = "neglog",
  alphaLink = "log",
  eimStep = 10
)

#y <- fn$simulate(n = N, eta = eta, lower = -1)
y <- testTruncated$ztnegbin3

df <- data.frame(
  x3 = x3[y > 0],
  x2 = x2[y > 0],
  x1 = x1[y > 0],
  y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = fn,
  data    = df,
  method = "IRLS",
  controlModel = controlModel(
    alphaFormula = ~ x1 + x2 + x3
  )
)

expect_silent(
  pop <- popSizeEst(M1)
)

expect_equivalent(
  N,
  pop$pointEstimate,
  tol = .6
)

expect_true(
  (N < pop$confidenceInterval$upperBound[2]) &
  (N > pop$confidenceInterval$lowerBound[2])
)

expect_silent(
  stra <- stratifyPopsize(M1, ~ x1)
)

expect_equivalent(
  as.numeric(table(df$x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Observed,
  tolerance = eps
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Estimated,
  tolerance = .6
)

expect_true(
  all(predict(
    M1, 
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    M1, 
    type = "link",
    se.fit = TRUE
  )[, 3:4]> 0)
)

expect_true(
  all(predict(
    M1, 
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)] >= 0)
)

expect_true(
  sum(predict(
    M1, 
    type = "contr"
  )) == M1$populationSize$pointEstimate
)

expect_silent(
  predict(
    M1, 
    type = "popSize",
    se.fit = TRUE
  )
)

# different link
beta <- c(.6, -.2, -1.25, -.1,
          .2, -.1, -.5,   .3)
#eta <- cbind(beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x3,
#             beta[5] + beta[6] * x1 + beta[7] * x2 + beta[8] * x3)

fn <- ztnegbin(
  lambdaLink = "log",
  alphaLink = "log",
  eimStep = 25
)

#y <- fn$simulate(n = N, eta = eta, lower = -1)
y <- testTruncated$ztnegbin4

df <- data.frame(
  x3 = x3[y > 0],
  x2 = x2[y > 0],
  x1 = x1[y > 0],
  y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = y ~ x1 + x2 + x3,
  model   = fn,
  data    = df,
  method = "IRLS",
  controlModel = controlModel(
    alphaFormula = ~ x1 + x2 + x3
  )
)

expect_equivalent(
  beta,
  coef(M1),
  tolerance = .4
)

expect_silent(
  pop <- popSizeEst(M1)
)

expect_equivalent(
  N,
  pop$pointEstimate,
  tol = .15
)

expect_true(
  all(N < pop$confidenceInterval$upperBound) &
  all(N > pop$confidenceInterval$lowerBound)
)

expect_silent(
  stra <- stratifyPopsize(M1, ~ x1)
)

expect_equivalent(
  as.numeric(table(df$x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Observed,
  tolerance = eps
)

expect_equivalent(
  as.numeric(table(x1))[as.numeric(substr(stra$name, start = 5, stop = 5)) + 1],
  stra$Estimated,
  tolerance = .15
)

AA <- summary(marginalFreq(M1), dropl5 = "group")

expect_true(
  all(AA$Test$`P(>X^2)` > .05)
)

expect_true(
  all(predict(
    M1, 
    type = "response",
    se.fit = TRUE
  )[, 3:4] > 0)
)

expect_true(
  all(predict(
    M1, 
    type = "link",
    se.fit = TRUE
  )[, 3:4]> 0)
)

expect_true(
  all(predict(
    M1, 
    type = "mean",
    se.fit = TRUE
  )[, c(2, 4)] >= 0)
)

expect_true(
  sum(predict(
    M1, 
    type = "contr"
  )) == M1$populationSize$pointEstimate
)

expect_silent(
  predict(
    M1, 
    type = "popSize",
    se.fit = TRUE
  )
)
