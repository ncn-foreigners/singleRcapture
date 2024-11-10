
eps <- if (capabilities("long.double")) sqrt(.Machine$double.eps) else 0.01

# test simulate
set.seed(123)
expect_equivalent(
  {N <- 10000
   gender <- rbinom(N, 1, 0.2)
   eta <- -1 + 0.5*gender
   counts <- rnbinom(N, mu = exp(eta), size = 1)
   df <- data.frame(gender, eta, counts)
   df2 <- subset(df, counts > 0)
   mod1 <-  estimatePopsize(
     formula = counts ~ 1 + gender,
     data = df2,
     model = "ztnegbin",
     method = "optim"
   )
   mid1_sim <- simulate(mod1, 10)
   dim(mid1_sim)
  },
  c(2920, 10),
  tolerance = eps
)

expect_silent(
  Model <- estimatePopsize(
    formula = capture ~ nation + age + gender, 
    data = netherlandsimmigrant, 
    model = ztpoisson, 
    method = "IRLS",
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_silent(
  Model1 <- estimatePopsize(
    formula = capture ~ 1, 
    data = netherlandsimmigrant, 
    model = ztpoisson, 
    method = "IRLS",
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_silent(
  Model2 <- estimatePopsize(
    formula = capture ~ . - age - reason, 
    data = netherlandsimmigrant, 
    model = zelterman, 
    method = "IRLS",
    controlMethod = controlMethod(silent = TRUE)
  )
)

df <- netherlandsimmigrant[, c(1:3,5)]
df$ww <- 0
### this is dplyr::count but slower and without dependencies
df <- aggregate(ww ~ ., df, FUN = length)

expect_silent(
  Model6 <- estimatePopsize(
    formula = capture ~ nation + age + gender, 
    data = df, 
    model = ztpoisson, 
    method = "IRLS",
    weights = df$ww,
    controlMethod = controlMethod(silent = TRUE),
    controlModel = controlModel(weightsAsCounts = TRUE)
  )
)

expect_equal(
  nobs(Model6),
  nobs(Model),
  tolerance = eps
)

expect_equal(
  Model$populationSize$pointEstimate,
  Model6$populationSize$pointEstimate,
  tolerance = eps
)

expect_equal(
  Model$populationSize$confidenceInterval,
  Model6$populationSize$confidenceInterval,
  tolerance = eps
)

expect_equal(
  Model$populationSize$variance,
  Model6$populationSize$variance,
  tolerance = eps
)

expect_equal(
  Model$coefficients,
  Model6$coefficients,
  tolerance = eps
)

expect_equal(
  Model$logL,
  Model6$logL,
  tolerance = eps
)

# dfbetas and dfpopsize
# 4 takes too long
expect_silent(
  dfb <- dfbeta(Model)
)

expect_silent(
  dfb1 <- dfbeta(Model1)
)

expect_silent(
  dfb2 <- dfbeta(Model2)
)

expect_silent(
  dfb6 <- dfbeta(Model6)
)

expect_silent(
  hatvalues(Model)
)

expect_silent(
  hatvalues(Model2)
)

expect_silent(
  dfp <- dfpopsize(Model, dfbeta = dfb)
)

expect_silent(
  dfp6 <- dfpopsize(Model6, dfbeta = dfb6)
)

expect_equal(
  max(abs(dfp)),
  4236.412,
  tolerance = .05
)

expect_equal(
  abs(mean(abs(dfp))),
  19.19,
  tolerance = .05
)

expect_true(
  abs(max(abs(dfpopsize(Model1, dfbeta = dfb1))) - 88.349) < .2
)

expect_true(
  abs(mean(abs(dfpopsize(Model1, dfbeta = dfb1))) - 8.1945) < .1
)

expect_true(
  abs(max(abs(dfpopsize(Model2, dfbeta = dfb2))) - 3648.17) < .2
)

expect_equal(
  c(unique(dfp[netherlandsimmigrant$capture == 1 & netherlandsimmigrant$gender == "female" & netherlandsimmigrant$nation == "American and Australia" & netherlandsimmigrant$age == "<40yrs"])[1],
    unique(dfp[netherlandsimmigrant$capture == 2 & netherlandsimmigrant$gender == "female" & netherlandsimmigrant$nation == "American and Australia" & netherlandsimmigrant$age == "<40yrs"])[1],
    unique(dfp[netherlandsimmigrant$capture == 3 & netherlandsimmigrant$gender == "female" & netherlandsimmigrant$nation == "American and Australia" & netherlandsimmigrant$age == "<40yrs"])[1],
    unique(dfp[netherlandsimmigrant$capture == 1 & netherlandsimmigrant$gender == "male"   & netherlandsimmigrant$nation == "American and Australia" & netherlandsimmigrant$age == "<40yrs"])[1],
    unique(dfp[netherlandsimmigrant$capture == 2 & netherlandsimmigrant$gender == "male"   & netherlandsimmigrant$nation == "American and Australia" & netherlandsimmigrant$age == "<40yrs"])[1],
    unique(dfp[netherlandsimmigrant$capture == 3 & netherlandsimmigrant$gender == "male"   & netherlandsimmigrant$nation == "American and Australia" & netherlandsimmigrant$age == "<40yrs"])[1],
    unique(dfp[netherlandsimmigrant$capture == 4 & netherlandsimmigrant$gender == "male"   & netherlandsimmigrant$nation == "American and Australia" & netherlandsimmigrant$age == "<40yrs"])[1]),
  dfp6[1:7],
  tolerance = eps
)

# Extractors

expect_true(
  max(abs(
    c(AIC(Model), AIC(Model1), AIC(Model2)) -
      c(1712.901, 1805.904, 1131.723)
  )) < 1
)

expect_true(
  max(abs(
    c(BIC(Model), BIC(Model1), BIC(Model2)) -
      c(1757.213, 1811.443, 1170.3)
  )) < 1
)

expect_silent(
  c(extractAIC(Model), extractAIC(Model1), extractAIC(Model2))
)

expect_equal(AIC(Model), AIC(Model6),
             tolerance = eps)
expect_equal(BIC(Model), BIC(Model6),
             tolerance = eps)
expect_equal(extractAIC(Model), extractAIC(Model6),
             tolerance = eps)

# Sandwich

require(sandwich)

expect_equivalent(
  vcovHC(Model),
  vcovHC(Model6),
  tolerance = eps
)

expect_equivalent(
  vcov(Model, type = "observedInform"),
  vcov(Model, type = "Fisher"),
  tolerance = max(.0001, eps)
)

expect_equivalent(
  vcov(Model1, type = "observedInform"),
  vcov(Model1, type = "Fisher"),
  tolerance = max(.0001, eps)
)

expect_equivalent(
  vcov(Model2, type = "observedInform"),
  vcov(Model2, type = "Fisher"),
  tolerance = max(.00001, eps)
)

expect_equivalent(
  bread(Model, type = "observedInform"),
  bread(Model, type = "Fisher"),
  tolerance = max(.0001, eps)
)

expect_equivalent(
  bread(Model1, type = "observedInform"),
  bread(Model1, type = "Fisher"),
  tolerance = max(.0001, eps)
)

expect_equivalent(
  bread(Model2, type = "observedInform"),
  bread(Model2, type = "Fisher"),
  tolerance = max(.00001, eps)
)

expect_silent(
  sandwich(Model)
)

expect_silent(
  sandwich(Model1)
)

expect_silent(
  sandwich(Model2)
)

expect_silent(
  sandwich(Model6)
)

expect_silent(
  vcovHC(Model)
)

expect_silent(
  vcovHC(Model6, type = "HC")
)

expect_silent(
  vcovHC(Model1, type = "HC4m")
)

expect_silent(
  vcovHC(Model2, type = "HC5")
)

expect_silent(
  vcovHC(Model, type = "const")
)

expect_silent(
  vcovHC(Model6, type = "HC2")
)

expect_silent(
  vcovHC(Model1, type = "HC1")
)

# confint

expect_silent(
  confint(Model)
)

expect_silent(
  confint(Model, parm = 1:2, level = .99)
)

expect_silent(
  cooks.distance(Model)
)

expect_silent(
  cooks.distance(Model1)
)

expect_silent(
  cooks.distance(Model2)
)

expect_silent(
  model.frame(Model2)
)

expect_true(
  all(dim(model.matrix(Model)) == c(1880, 8))
)

temp <- Model
temp$X <- NULL

expect_identical(
  model.matrix(Model),
  model.matrix(temp)
)

expect_identical(
  model.matrix(Model, "vlm"),
  model.matrix(temp, "vlm")
)

temp$modelFrame <- NULL

expect_identical(
  model.matrix(Model),
  model.matrix(temp)
)

expect_identical(
  model.matrix(Model, "vlm"),
  model.matrix(temp, "vlm")
)

expect_silent(
  print(popSizeEst(Model))
)

expect_equal(
  c(popSizeEst(Model)$pointEstimate, popSizeEst(Model1)$pointEstimate, 
    popSizeEst(Model2)$pointEstimate),
  c(12690, 7080, 15816.14),
  tol = .05
)

expect_equal(
  c(popSizeEst(Model)$variance, popSizeEst(Model1)$variance, 
    popSizeEst(Model2)$variance),
  c(7885812, 133774.1, 9093464),
  tol = .05
)

expect_silent(
  plot(Model, "qq")
)

expect_silent(
  plot(Model, "marginal")
)

expect_silent(
  plot(Model, "fitresid")
)

expect_error(
  plot(Model, "bootHist")
)

expect_silent(
  plot(Model, "rootogram")
)

expect_silent(
  plot(Model, "dfpopContr", dfpop = dfp)
)

expect_silent(
  plot(Model, "dfpopBox", dfpop = dfp)
)

expect_silent(
  plot(Model, "scaleLoc")
)

expect_silent(
  plot(Model, "cooks")
)

expect_silent(
  plot(Model, "hatplot")
)

expect_silent(
  plot(Model, "qq")
)

expect_silent(
  plot(Model, "strata")
)

expect_silent(
  up <- redoPopEstimation(Model, cov = vcovHC(Model, "HC4m"))
)

expect_silent(
  summary(Model, cov = vcovHC(Model, "HC4m"), correlation = TRUE, 
          confint = TRUE, popSizeEst = up)
)

expect_equivalent(
  up$confidenceInterval[1, ],
  data.frame(6611.906, 18768.8),
  tolerance = .05
)

expect_silent(
  up <- summary(marginalFreq(Model6), df = 1, dropl5 = "group")
)

expect_silent(
  print(up)
)

expect_equivalent(
  predict(Model, type = "response", se.fit = TRUE),
  predict(Model6, type = "response", se.fit = TRUE, newdata = model.frame(Model)),
  tolerance = eps
)

expect_equivalent(
  predict(Model, type = "link", se.fit = TRUE),
  predict(Model6, type = "link", se.fit = TRUE, newdata = netherlandsimmigrant[,-4]),
  tolerance = eps
)

expect_equivalent(
  predict(Model, type = "mean", se.fit = TRUE),
  predict(Model6, type = "mean", se.fit = TRUE, newdata = netherlandsimmigrant[,-4]),
  tolerance = eps
)

expect_silent(
  residuals(Model, type = "all")
)

expect_equal(
  logLik(Model),
  logLik(Model6),
  tolerance = eps
)

expect_silent(
  up <- redoPopEstimation(
    Model6, newdata = netherlandsimmigrant[,-4]
  )
)

expect_silent(
  up1 <- redoPopEstimation(
    Model, 
    newdata = model.frame(Model6), 
    weights = Model6$priorWeights,
    weightsAsCounts = TRUE
  )
)

expect_equal(
  up, up1, tolerance = .025
)

expect_equal(
  stratifyPopsize(Model),
  stratifyPopsize(Model6),
  tolerance = eps
)

expect_error(
  stratifyPopsize(Model, strata = 8L)
)

expect_silent(
  print(Model)
)

expect_silent(
  print(family(Model))
)

expect_equal(
  NCOL(fitted(Model, "all")),
  2L
)

expect_equivalent(
  as.numeric(table(simulate(Model6, seed = 123)[,1])),
  c(1619, 232, 28, 1),
  tolerance = eps
)
