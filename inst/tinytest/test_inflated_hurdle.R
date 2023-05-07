#### tracking progress:
##
## - ztoigeom - test
## - oiztgeom - test
## - ztHurdlegeom - test
## - Hurdleztgeom - test
## - ztoipoisson - test
## - oiztpoisson - test
## - ztHurdlepoisson - test
## - Hurdleztpoisson - test
## - ztoinegbin - to be finished -- variance, W, test
## - oiztnegbin - to started
## - ztHurdlenegbin - to be started
## - Hurdleztnegbin - to be started

set.seed(123)

N <- 10000
x1 <- rbinom(n = N, size = 3, prob = .4)


# Ztoi ####

# ztoipoisson
beta <- c(.6, .2, -1.25, .1)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)

fn <- ztoipoisson(
  lambdaLink = "log",
  omegaLink  = "logit",
)

y <- fn$simulate(n = N, eta = eta, lower = -1)

df <- data.frame(
  X = x1[y > 0],
  Y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = Y ~ X,
  model   = fn,
  data    = df,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ X
  )
)

expect_true(
  all(abs(coef(M1) - beta) < .5)
)

expect_true(
  abs(popSizeEst(M1)$pointEstimate - N) < 100
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
  (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
  (N < popSizeEst(M1)$confidenceInterval[2, 2])
)

# same with different link

fn <- ztoipoisson(
  lambdaLink = "neglog",
  omegaLink  = "cloglog",
)

y <- fn$simulate(n = N, eta = eta, lower = -1)

df <- data.frame(
  X = x1[y > 0],
  Y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = Y ~ X,
  model   = fn,
  data    = df,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ X
  )
)

expect_true(
  all(abs(coef(M1) - beta) < 1.15)
)

expect_true(
  abs(popSizeEst(M1)$pointEstimate - N) < 1500
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[2, 2])
)

# ztoigeom
beta <- c(.6, .2, -1.25, .1)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)

fn <- ztoigeom(
  lambdaLink = "log",
  omegaLink  = "logit",
)

y <- fn$simulate(n = N, eta = eta, lower = -1)

df <- data.frame(
  X = x1[y > 0],
  Y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = Y ~ X,
  model   = fn,
  data    = df,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ X
  )
)

expect_true(
  all(abs(coef(M1) - beta) < .5)
)

expect_true(
  abs(popSizeEst(M1)$pointEstimate - N) < 180
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
  (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
  (N < popSizeEst(M1)$confidenceInterval[2, 2])
)


# same for different link

fn <- ztoigeom(
  lambdaLink = "neglog",
  omegaLink  = "cloglog",
)

y <- fn$simulate(n = N, eta = eta, lower = -1)

df <- data.frame(
  X = x1[y > 0],
  Y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = Y ~ X,
  model   = fn,
  data    = df,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ X
  )
)

expect_true(
  all(abs(coef(M1) - beta) < .5)
)

expect_true(
  abs(popSizeEst(M1)$pointEstimate - N) < 1000
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[2, 2])
)


# ztoinegbin
beta <- c(.6, .2, -1.25, .1, -1, -.3)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1, beta[5] + beta[6] * x1)

fn <- ztoinegbin(
  lambdaLink = "log",
  alphaLink  = "log", 
  omegaLink  = "logit",
)

y <- fn$simulate(n = N, eta = eta, lower = -1)

df <- data.frame(
  X = x1[y > 0],
  Y = y[y > 0]
)

M1 <- estimatePopsize(
  formula = Y ~ X,
  model   = fn,
  data    = df,
  method = "IRLS",
  controlModel = controlModel(
    omegaFormula = ~ X,
    alphaFormula = ~ X
  )
)

expect_true(
  all(abs(coef(M1) - beta) < .5)
)

expect_true(
  abs(popSizeEst(M1)$pointEstimate - N) < 400
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
  (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
  (N < popSizeEst(M1)$confidenceInterval[2, 2])
)

# oizt ####
