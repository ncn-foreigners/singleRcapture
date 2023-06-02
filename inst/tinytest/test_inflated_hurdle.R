# Testing on generated data for inflated and hurdle models

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

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

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

# test whether estimation of subpopulation sizes are 'good'
expect_equivalent(
  stratifyPopsize(M1, stratas = ~ X)$Estimated[c(4, 1:3)],
  as.numeric(table(x1)),
  tolerance = .15
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

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

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

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

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

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

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

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

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

# poisson
beta <- c(.6, -.3, -.5, .1)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)

fn <- oiztpoisson(
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

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

expect_true(
  all(abs(coef(M1) - beta) < .5)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .1
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

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

expect_true(
  all(abs(coef(M1) - beta) < .6)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .15
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
  (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[2, 2])
)

# geometric
beta <- c(.6, .2, -1, .4)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)

fn <- oiztgeom(
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

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

expect_true(
  all(abs(coef(M1) - beta) < .5)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .1
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

beta <- c(.6, -.3, -1, .4)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)


fn <- oiztgeom(
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

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

expect_true(
  all(abs(coef(M1) - beta) < .5)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .1
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[2, 2])
)


# negbin
beta <- c(.4, .4, -.4, .2, -1.1, .37)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1, beta[5] + beta[6] * x1)

fn <- oiztnegbin(
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

expect_silent(summary(M1, confint = TRUE, correlation = TRUE))

expect_silent(print(summary(M1, confint = TRUE, correlation = TRUE)))

expect_true(
  all(abs(coef(M1) - beta) < .6)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .1
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[2, 2])
)

# ztHurdle ####


# poisson
beta <- c(.6, -.3, -.5, .1)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)

fn <- ztHurdlepoisson(
  lambdaLink = "log",
  piLink  = "logit",
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
    piFormula = ~ X
  )
)

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

expect_true(
  all(abs(coef(M1) - beta) < .5)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .1
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

fn <- ztHurdlepoisson(
  lambdaLink = "neglog",
  piLink  = "probit",
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
    piFormula = ~ X
  )
)

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

expect_true(
  all(abs(coef(M1) - beta) < .3)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .15
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[2, 2])
)

# geometric
beta <- c(.6, .2, -3, .4)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)

fn <- ztHurdlegeom(
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
    piFormula = ~ X
  )
)

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

expect_true(
  all(abs(coef(M1) - beta) < .4)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .1
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

beta <- c(.6, -.3, -1, .4)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)


fn <- ztHurdlegeom(
  lambdaLink = "neglog", 
  piLink = "probit",
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
    piFormula = ~ X
  )
)

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

expect_true(
  all(abs(coef(M1) - beta) < .4)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .1
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[2, 2])
)


# negbin
beta <- c(.4, .4, -.4, .2, -1.1, .37)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1, beta[5] + beta[6] * x1)

fn <- ztHurdlenegbin(
  lambdaLink = "log",
  alphaLink  = "log", 
  piLink  = "logit",
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
    piFormula = ~ X,
    alphaFormula = ~ X
  )
)

expect_silent(summary(M1, confint = TRUE, correlation = TRUE))

expect_silent(print(summary(M1, confint = TRUE, correlation = TRUE)))

expect_true(
  all(abs(coef(M1) - beta) < .6)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .1
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[2, 2])
)
# Hurdlezt ####


# poisson
beta <- c(.6, -.3, -.5, .1)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)

fn <- Hurdleztpoisson(
  lambdaLink = "log",
  piLink  = "logit",
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
    piFormula = ~ X
  )
)

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

expect_true(
  all(abs(coef(M1) - beta) < .5)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .1
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

fn <- Hurdleztpoisson(
  lambdaLink = "neglog",
  piLink  = "probit",
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
    piFormula = ~ X
  )
)

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

expect_true(
  all(abs(coef(M1) - beta) < .3)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .15
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[2, 2])
)

# geometric
beta <- c(.6, .2, -3, .4)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)

fn <- Hurdleztgeom(
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
    piFormula = ~ X
  )
)

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

expect_true(
  all(abs(coef(M1) - beta) < .4)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .1
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

beta <- c(.3, -.7, -.25, .3)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)


fn <- Hurdleztgeom(
  lambdaLink = "neglog", 
  piLink = "probit",
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
    piFormula = ~ X
  )
)

expect_silent(summary(M1))

expect_silent(print(summary(M1)))

expect_true(
  all(abs(coef(M1) - beta) < .4)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .1
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[2, 2])
)


# negbin
beta <- c(.4, .4, -.4, .2, -1.1, .37)
eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1, beta[5] + beta[6] * x1)

fn <- Hurdleztnegbin(
  lambdaLink = "log",
  alphaLink  = "log", 
  piLink  = "logit",
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
    piFormula = ~ X,
    alphaFormula = ~ X
  )
)

expect_silent(summary(M1, confint = TRUE, correlation = TRUE))

expect_silent(print(summary(M1, confint = TRUE, correlation = TRUE)))

expect_true(
  all(abs(coef(M1) - beta) < .6)
)

expect_equal(
  popSizeEst(M1)$pointEstimate,
  N, tolerance = .1
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[1, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[1, 2])
)

expect_true(
  (popSizeEst(M1)$confidenceInterval[2, 1] < N) &
    (N < popSizeEst(M1)$confidenceInterval[2, 2])
)
