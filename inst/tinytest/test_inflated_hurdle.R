# Testing on generated data for inflated and hurdle models
# These tests are only supposed to be run on developer's machine and 
# package GitHub page not on CRAN (they take too long)

if (isTRUE(tolower(Sys.getenv("TEST_SINGLERCAPTURE_MULTICORE_DEVELOPER")) == "true")) {
    
  N <- 500
  #x1 <- rbinom(n = N, size = 3, prob = .4)
  x1 <- rep(0:3, c(99, 228, 137, 36))
  ## generated data from distributions
  test_inflated <- read.csv("test_inflated.csv")
  #test_inflated <- read.csv("inst/tinytest/test_inflated.csv")
  
  # Ztoi ####
  
  # ztoipoisson
  beta <- c(.6, .2, -1.25, .1)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)
  
  fn <- ztoipoisson(
    lambdaLink = "log",
    omegaLink  = "logit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$ztoipoisson1
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      omegaFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1), 
    beta, tol = .3
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate, 
    N, tol = .1
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
  )
  
  # test whether estimation of subpopulation sizes are 'good'
  XX <- stratifyPopsize(M1, stratas = ~ x)
  order <- 1 + as.numeric(substr(XX$name, start = 4, stop = 4))
  
  expect_equivalent(
    XX$Estimated,
    as.numeric(table(x1))[order],
    tolerance = .15
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # same with different link
  fn <- ztoipoisson(
    lambdaLink = "neglog",
    omegaLink  = "cloglog",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$ztoipoisson2
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      omegaFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .7
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate, N,
    tolerance = .25
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # ztoigeom
  beta <- c(.6, .2, -1.25, .1)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)
  
  fn <- ztoigeom(
    lambdaLink = "log",
    omegaLink  = "logit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$ztoigeom1
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      omegaFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .7
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .2
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # same for different link
  fn <- ztoigeom(
    lambdaLink = "neglog",
    omegaLink  = "cloglog",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$ztoigeom2
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      omegaFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .2
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # ztoinegbin
  beta <- c(.6, .2, -1.25, .1, -1, -.3)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1, beta[5] + beta[6] * x1)
  
  fn <- ztoinegbin(
    lambdaLink = "log",
    alphaLink  = "log", 
    omegaLink  = "logit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$ztoinegbin1
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      omegaFormula = ~ x,
      alphaFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  # don't change this
  expect_true(
    sum(abs(coef(M1) - beta)) < 8
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .3
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "response",
      se.fit = TRUE
    )[, 4:6] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "link",
      se.fit = TRUE
    )[, 4:6]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 4:6] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[, 4:6]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # oizt ####
  
  # poisson
  beta <- c(.6, -.3, -.5, .1)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)
  
  fn <- oiztpoisson(
    lambdaLink = "log",
    omegaLink  = "logit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$oiztpoisson1
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      omegaFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .25
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .2
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  
  # same with different link
  fn <- oiztpoisson(
    lambdaLink = "neglog",
    omegaLink  = "cloglog",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$oiztpoisson2
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      omegaFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .5
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .2
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # geometric
  beta <- c(.6, .2, -1, .4)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)
  
  fn <- oiztgeom(
    lambdaLink = "log",
    omegaLink  = "logit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$oiztgeom1
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      omegaFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .3
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .2
  )
  
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # same for different link
  beta <- c(.6, -.3, -1, .4)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)
  
  fn <- oiztgeom(
    lambdaLink = "neglog",
    omegaLink  = "cloglog",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$oiztgeom2
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      omegaFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .7
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .2
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # negbin
  beta <- c(.4, .4, -.4, .2, -1.1, .37)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1, beta[5] + beta[6] * x1)
  
  fn <- oiztnegbin(
    lambdaLink = "log",
    alphaLink  = "log", 
    omegaLink  = "logit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$oiztnegbin1
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      omegaFormula = ~ x,
      alphaFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1, confint = TRUE, correlation = TRUE))
  
  expect_silent(print(summary(M1, confint = TRUE, correlation = TRUE)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .6
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .2
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "response",
      se.fit = TRUE
    )[, 4:6] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "link",
      se.fit = TRUE
    )[, 4:6]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 4:6] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[, 4:6]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # ztHurdle ####
  
  # poisson
  beta <- c(.6, -.3, -.5, .1)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)
  
  fn <- ztHurdlepoisson(
    lambdaLink = "log",
    piLink  = "logit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$ztHurdlepoisson1
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      piFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .7
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .5
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # same with different link
  fn <- ztHurdlepoisson(
    lambdaLink = "neglog",
    piLink  = "probit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$ztHurdlepoisson2
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      piFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .4
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .4
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # geometric
  beta <- c(.6, .2, -3, .4)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)
  
  fn <- ztHurdlegeom(
    lambdaLink = "log",
    omegaLink  = "logit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$ztHurdlegeom1
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      piFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .4
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .2
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # same for different link
  beta <- c(.6, -.3, -1, .4)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)
  
  
  fn <- ztHurdlegeom(
    lambdaLink = "neglog", 
    piLink = "probit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$ztHurdlegeom2
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      piFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .65
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .25
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # negbin
  beta <- c(.4, .4, -.4, .2, -1.1, .37)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1, beta[5] + beta[6] * x1)
  
  fn <- ztHurdlenegbin(
    lambdaLink = "log",
    alphaLink  = "log", 
    piLink  = "logit",
    eimStep = 20
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$ztHurdlenegbin1
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      piFormula = ~ x,
      alphaFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1, confint = TRUE, correlation = TRUE))
  
  expect_silent(print(summary(M1, confint = TRUE, correlation = TRUE)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .5
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .3
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "response",
      se.fit = TRUE
    )[, 4:6] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "link",
      se.fit = TRUE
    )[, 4:6]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 4:6] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[, 4:6]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # Hurdlezt ####
  # poisson
  beta <- c(.6, -.3, -.5, .1)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)
  
  fn <- Hurdleztpoisson(
    lambdaLink = "log",
    piLink  = "logit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$Hurdleztpoisson1
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      piFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .6
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .25
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # same with different link
  fn <- Hurdleztpoisson(
    lambdaLink = "neglog",
    piLink  = "probit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$Hurdleztpoisson2
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      piFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE, stepsize = .3)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .6
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .3
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # geometric
  beta <- c(.6, .2, -3, .4)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)
  
  fn <- Hurdleztgeom(
    lambdaLink = "log",
    omegaLink  = "logit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$Hurdleztgeom1
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      piFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .2
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .2
  )
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
  )
  
  # same for different link
  beta <- c(.3, -.7, -.25, .3)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1)
  
  fn <- Hurdleztgeom(
    lambdaLink = "neglog", 
    piLink = "probit",
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$Hurdleztgeom2
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      piFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1))
  
  expect_silent(print(summary(M1)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .4
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .2
  )
  
  
  expect_true(
    (pop$confidenceInterval[1, 1] < N) &
    (N < pop$confidenceInterval[1, 2])
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
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
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 3:4] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[,3:4]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  # negbin
  beta <- c(.4, .4, -.4, .2, -1.1, .37)
  #eta <- cbind(beta[1] + beta[2] * x1, beta[3] + beta[4] * x1, beta[5] + beta[6] * x1)
  
  fn <- Hurdleztnegbin(
    lambdaLink = "log",
    alphaLink  = "log", 
    piLink  = "logit",
    eimStep = 30
  )
  
  #y <- fn$simulate(n = N, eta = eta, lower = -1)
  y <- test_inflated$Hurdleztnegbin1
  
  df <- data.frame(
    x = x1[y > 0],
    y = y[y > 0]
  )
  
  M1 <- estimatePopsize(
    formula = y ~ x,
    model   = fn,
    data    = df,
    method = "IRLS",
    controlModel = controlModel(
      piFormula = ~ x,
      alphaFormula = ~ x
    ),
    controlMethod = controlMethod(silent = TRUE)
  )
  
  expect_silent(summary(M1, confint = TRUE, correlation = TRUE))
  
  expect_silent(print(summary(M1, confint = TRUE, correlation = TRUE)))
  
  expect_equivalent(
    coef(M1),
    beta,
    tolerance = .7
  )
  
  expect_silent(
    pop <- popSizeEst(M1)
  )
  
  expect_equivalent(
    pop$pointEstimate,
    N,
    tolerance = .35
  )
  
  expect_true(
    (pop$confidenceInterval[2, 1] < N) &
    (N < pop$confidenceInterval[2, 2])
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "response",
      se.fit = TRUE
    )[, 4:6] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "link",
      se.fit = TRUE
    )[, 4:6]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_equivalent(
    sum(predict(
      M1, 
      type = "contr"
    )),
    M1$populationSize$pointEstimate
  )
  
  expect_silent(
    predict(
      M1, 
      type = "popSize",
      se.fit = TRUE
    )
  )
  
  expect_true(
    all(predict(
      M1,
      newdata = data.frame(x = x1),
      type = "response",
      se.fit = TRUE
    )[, 4:6] > 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "link",
      se.fit = TRUE
    )[, 4:6]> 0)
  )
  
  expect_true(
    all(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "mean",
      se.fit = TRUE
    )[, c(2, 4)]> 0)
  )
  
  expect_true(
    sum(predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "contr"
    )) > N
  )
  
  expect_silent(
    predict(
      M1, 
      newdata = data.frame(x = x1),
      type = "popSize",
      se.fit = TRUE
    )
  )
}