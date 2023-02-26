# test simulate
# set.seed(123)
# expect_equivalent(
#   {N <- 10000
#    gender <- rbinom(N, 1, 0.2)
#    eta <- -1 + 0.5*gender
#    disp <- 1
#    counts <- rnbinom(N, mu = exp(eta), size = disp)
#    df <- data.frame(gender, eta, counts)
#    df2 <- subset(df, counts > 0)
#    mod1 <-  estimatePopsize(formula = counts ~ 1 + gender, 
#                              data = df2,
#                              model = "ztnegbin", 
#                              method = "optim",
#                              pop.var = "analytic")
#    mid1_sim <- simulate(mod1, 10)
#    dim(mid1_sim)
#   },
#   c(2920, 10)
# )

# we only check methods in this file
Model <- estimatePopsize(
  formula = capture ~ nation + age + gender, 
  data = netherlandsimmigrant, 
  model = ztpoisson, 
  method = "IRLS"
)

Model1 <- estimatePopsize(
  formula = capture ~ 1, 
  data = netherlandsimmigrant, 
  model = ztpoisson, 
  method = "IRLS"
)

Model2 <- estimatePopsize(
  formula = capture ~ . - age - reason, 
  data = netherlandsimmigrant, 
  model = zelterman, 
  method = "IRLS"
)

Model3 <- estimatePopsize(
  formula = capture ~ . - age, 
  data = netherlandsimmigrant, 
  model = chao, 
  method = "IRLS"
)

Model4 <- estimatePopsize(
  formula = TOTAL_SUB ~ ., 
  data = farmsubmission, 
  model = ztnegbin,
  method = "IRLS"
)

Model5 <- estimatePopsize(
  formula = TOTAL_SUB ~ ., 
  data = farmsubmission, 
  model = ztoigeom, 
  method = "IRLS",
  controlPopVar = controlPopVar(covType = "Fisher"),
  controlModel = controlModel(omegaFormula = ~ log_distance),
  controlMethod = controlMethod(stepsize = .45)
)

# dfbetas and dfpopsize

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
  dfb3 <- dfbeta(Model3)
)

expect_equivalent(
  max(abs(dfpopsize(Model, dfbeta = dfb))),
  4236.412,
  tolerance = .1
)

expect_equivalent(
  max(abs(dfpopsize(Model1, dfbeta = dfb1))),
  88.35,
  tolerance = .1
)

expect_equivalent(
  max(abs(dfpopsize(Model2, dfbeta = dfb2))),
  3648.162,
  tolerance = .1
)

expect_equivalent(
  max(abs(dfpopsize(Model3, dfbeta = dfb3))),
  3681.764,
  tolerance = .1
)

# Extractors

expect_equivalent(
  AIC(Model),
  1712.901,
  tolerance = .1
)

expect_equivalent(
  AIC(Model1),
  1805.904,
  tolerance = .1
)

expect_equivalent(
  AIC(Model2),
  1131.723,
  tolerance = .1
)

expect_equivalent(
  AIC(Model3),
  1131.006,
  tolerance = .1
)

expect_equivalent(
  AIC(Model4),
  34538.73,
  tolerance = .1
)

expect_equivalent(
  AIC(Model5),
  34580.82,
  tolerance = .1
)

expect_equivalent(
  BIC(Model),
  1757.213,
  tolerance = .1
)

expect_equivalent(
  BIC(Model1),
  1811.443,
  tolerance = .1
)

expect_equivalent(
  BIC(Model2),
  1170.3,
  tolerance = .1
)

expect_equivalent(
  BIC(Model3),
  1177.094,
  tolerance = .1
)

expect_equivalent(
  BIC(Model4),
  34575.71,
  tolerance = .1
)

expect_equivalent(
  BIC(Model5),
  34625.2,
  tolerance = .1
)

expect_silent(
  extractAIC(Model)
)

expect_silent(
  extractAIC(Model1)
)

expect_silent(
  extractAIC(Model2)
)

expect_silent(
  extractAIC(Model3)
)

expect_silent(
  extractAIC(Model4)
)

expect_silent(
  extractAIC(Model5)
)

# Sandwich

require(sandwich)

expect_silent(
  bread(Model)
)

expect_silent(
  bread(Model1)
)

expect_silent(
  bread(Model2)
)

expect_silent(
  bread(Model3)
)

expect_silent(
  bread(Model4)
)

expect_silent(
  bread(Model5)
)

expect_silent(
  bread(Model, type = "Fisher")
)

expect_silent(
  bread(Model1, type = "Fisher")
)

expect_silent(
  bread(Model2, type = "Fisher")
)

expect_silent(
  bread(Model3, type = "Fisher")
)

expect_silent(
  bread(Model4, type = "Fisher")
)

expect_silent(
  bread(Model5, type = "Fisher")
)

expect_false(
  all(vcov(Model, type = "observedInform") == vcov(Model, type = "Fisher"))
)

expect_false(
  all(vcov(Model1, type = "observedInform") == vcov(Model1, type = "Fisher"))
)

expect_false(
  all(vcov(Model4, type = "observedInform") == vcov(Model4, type = "Fisher"))
)

expect_false(
  all(vcov(Model5, type = "observedInform") == vcov(Model5, type = "Fisher"))
)

expect_equal(
  sum(vcov(Model2, type = "observedInform")-vcov(Model2, type = "Fisher")),
  0, tolerance = 1e-5
)

expect_true(
  max(vcov(Model3, type = "Fisher") - vcov(Model3, type = "observedInform")) < 1e-8
)

expect_false(
  all(bread(Model, type = "observedInform") == bread(Model, type = "Fisher"))
)

expect_false(
  all(bread(Model1, type = "observedInform") == bread(Model1, type = "Fisher"))
)

expect_true(
  max(bread(Model2, type = "observedInform") - bread(Model2, type = "Fisher")) < 1e-1
)

expect_true(
  max(bread(Model3, type = "observedInform") - bread(Model3, type = "Fisher")) < 1e-1
)

expect_false(
  all(bread(Model4, type = "observedInform") == bread(Model4, type = "Fisher"))
)

expect_false(
  all(bread(Model5, type = "observedInform") == bread(Model5, type = "Fisher"))
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
  sandwich(Model3)
)

expect_silent(
  sandwich(Model4)
)

expect_silent(
  sandwich(Model5)
)

expect_silent(
  vcovHC(Model)
)

expect_silent(
  vcovHC(Model1)
)

expect_silent(
  vcovHC(Model2)
)

expect_silent(
  vcovHC(Model3)
)

expect_silent(
  vcovHC(Model4)
)

expect_silent(
  vcovHC(Model5)
)

expect_silent(
  vcovHC(Model, type = "HC")
)

expect_silent(
  vcovHC(Model1, type = "HC")
)

expect_silent(
  vcovHC(Model2, type = "HC")
)

expect_silent(
  vcovHC(Model3, type = "HC")
)

expect_silent(
  vcovHC(Model4, type = "HC")
)

expect_silent(
  vcovHC(Model5, type = "HC")
)

expect_silent(
  vcovHC(Model, type = "HC0")
)

expect_silent(
  vcovHC(Model1, type = "HC0")
)

expect_silent(
  vcovHC(Model2, type = "HC0")
)

expect_silent(
  vcovHC(Model3, type = "HC0")
)

expect_silent(
  vcovHC(Model4, type = "HC0")
)

expect_silent(
  vcovHC(Model5, type = "HC0")
)

expect_silent(
  vcovHC(Model, type = "HC1")
)

expect_silent(
  vcovHC(Model1, type = "HC1")
)

expect_silent(
  vcovHC(Model2, type = "HC1")
)

expect_silent(
  vcovHC(Model3, type = "HC1")
)

expect_silent(
  vcovHC(Model4, type = "HC1")
)

expect_silent(
  vcovHC(Model5, type = "HC1")
)

expect_silent(
  vcovHC(Model, type = "HC2")
)

expect_silent(
  vcovHC(Model1, type = "HC2")
)

expect_silent(
  vcovHC(Model2, type = "HC2")
)

expect_silent(
  vcovHC(Model3, type = "HC2")
)

expect_silent(
  vcovHC(Model4, type = "HC2")
)

expect_silent(
  vcovHC(Model5, type = "HC2")
)

expect_silent(
  vcovHC(Model, type = "HC3")
)

expect_silent(
  vcovHC(Model1, type = "HC3")
)

expect_silent(
  vcovHC(Model2, type = "HC3")
)

expect_silent(
  vcovHC(Model3, type = "HC3")
)

expect_silent(
  vcovHC(Model4, type = "HC3")
)

expect_silent(
  vcovHC(Model5, type = "HC3")
)

expect_silent(
  vcovHC(Model, type = "HC4")
)

expect_silent(
  vcovHC(Model1, type = "HC4")
)

expect_silent(
  vcovHC(Model2, type = "HC4")
)

expect_silent(
  vcovHC(Model3, type = "HC4")
)

expect_silent(
  vcovHC(Model4, type = "HC4")
)

expect_silent(
  vcovHC(Model5, type = "HC4")
)

expect_silent(
  vcovHC(Model, type = "HC4m")
)

expect_silent(
  vcovHC(Model1, type = "HC4m")
)

expect_silent(
  vcovHC(Model2, type = "HC4m")
)

expect_silent(
  vcovHC(Model3, type = "HC4m")
)

expect_silent(
  vcovHC(Model4, type = "HC4m")
)

expect_silent(
  vcovHC(Model5, type = "HC4m")
)

expect_silent(
  vcovHC(Model, type = "HC5")
)

expect_silent(
  vcovHC(Model1, type = "HC5")
)

expect_silent(
  vcovHC(Model2, type = "HC5")
)

expect_silent(
  vcovHC(Model3, type = "HC5")
)

expect_silent(
  vcovHC(Model4, type = "HC5")
)

expect_silent(
  vcovHC(Model5, type = "HC5")
)

expect_silent(
  vcovHC(Model, type = "const")
)

expect_silent(
  vcovHC(Model1, type = "const")
)

expect_silent(
  vcovHC(Model2, type = "const")
)

expect_silent(
  vcovHC(Model3, type = "const")
)

expect_silent(
  vcovHC(Model4, type = "const")
)

expect_silent(
  vcovHC(Model5, type = "const")
)

# confint

expect_silent(
  confint(Model)
)

expect_silent(
  confint(Model1)
)

expect_silent(
  confint(Model2)
)

expect_silent(
  confint(Model3)
)

expect_silent(
  confint(Model4)
)

expect_silent(
  confint(Model5)
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
  cooks.distance(Model3)
)

expect_error(
  cooks.distance(Model4)
)

expect_error(
  cooks.distance(Model5)
)

expect_identical(
  family(Model),
  Model$model
)

expect_identical(
  family(Model1),
  Model1$model
)

expect_identical(
  family(Model2),
  Model2$model
)

expect_identical(
  family(Model3),
  Model3$model
)

expect_silent(
  model.frame(Model)
)

expect_silent(
  model.frame(Model1)
)

expect_silent(
  model.frame(Model2)
)

expect_silent(
  model.frame(Model3)
)

expect_silent(
  model.frame(Model4)
)

expect_silent(
  model.frame(Model5)
)

expect_silent(
  model.matrix(Model)
)

expect_silent(
  model.matrix(Model1)
)

expect_silent(
  model.matrix(Model2)
)

expect_silent(
  model.matrix(Model3)
)

expect_silent(
  model.matrix(Model4)
)

expect_silent(
  model.matrix(Model5)
)

expect_true(
  all(dim(model.matrix(Model)) == c(1880, 8))
)

expect_true(
  all(dim(model.matrix(Model1)) == c(1880, 1))
)

expect_true(
  all(dim(model.matrix(Model2)) == c(1828, 7))
)

expect_true(
  all(dim(model.matrix(Model3)) == c(1828, 8))
)

expect_true(
  all(dim(model.matrix(Model4)) == c(12036, 4))
)

expect_true(
  all(dim(model.matrix(Model5)) == c(12036, 4))
)

expect_true(
  all(dim(model.matrix(Model4, "vlm")) == c(2 * 12036, 5))
)

expect_true(
  all(dim(model.matrix(Model5, "vlm")) == c(2 * 12036, 6))
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
  popSizeEst(Model)
)

expect_silent(
  popSizeEst(Model1)
)

expect_silent(
  popSizeEst(Model2)
)

expect_silent(
  popSizeEst(Model3)
)

expect_silent(
  popSizeEst(Model4)
)

expect_silent(
  popSizeEst(Model5)
)

expect_equal(
  popSizeEst(Model)$pointEstimate,
  12690,
  tolerance = .5
)

expect_equal(
  popSizeEst(Model1)$pointEstimate,
  7080,
  tolerance = .5
)

expect_equal(
  popSizeEst(Model2)$pointEstimate,
  15816.14,
  tolerance = .5
)

expect_equal(
  popSizeEst(Model3)$pointEstimate,
  15713.14,
  tolerance = .5
)

expect_equal(
  popSizeEst(Model4)$pointEstimate,
  41020.3,
  tolerance = .5
)

expect_equal(
  popSizeEst(Model5)$pointEstimate,
  29478.72,
  tolerance = .5
)

expect_equal(
  popSizeEst(Model)$variance,
  7885812,
  tolerance = .5
)

expect_equal(
  popSizeEst(Model1)$variance,
  133774.1,
  tolerance = .5
)

expect_equal(
  popSizeEst(Model2)$variance,
  9093464,
  tolerance = .5
)

expect_equal(
  popSizeEst(Model3)$variance,
  9096077,
  tolerance = .5
)

expect_equal(
  popSizeEst(Model4)$variance,
  3563463,
  tolerance = .5
)

expect_equal(
  popSizeEst(Model5)$variance,
  431426.9,
  tolerance = .5
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
  plot(Model, "dfpopContr")
)

expect_silent(
  plot(Model, "dfpopBox")
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
  plot(Model1, "qq")
)

expect_silent(
  plot(Model1, "marginal")
)

expect_silent(
  plot(Model1, "fitresid")
)

expect_error(
  plot(Model1, "bootHist")
)

expect_silent(
  plot(Model1, "rootogram")
)

expect_silent(
  plot(Model1, "dfpopContr")
)

expect_silent(
  plot(Model1, "dfpopBox")
)

expect_silent(
  plot(Model1, "scaleLoc")
)

expect_silent(
  plot(Model1, "cooks")
)

expect_silent(
  plot(Model1, "hatplot")
)

expect_silent(
  plot(Model2, "qq")
)

expect_silent(
  plot(Model2, "marginal")
)

expect_silent(
  plot(Model2, "fitresid")
)

expect_error(
  plot(Model2, "bootHist")
)

expect_silent(
  plot(Model2, "rootogram")
)

expect_silent(
  plot(Model2, "dfpopContr")
)

expect_silent(
  plot(Model2, "dfpopBox")
)

expect_silent(
  plot(Model2, "scaleLoc")
)

expect_silent(
  plot(Model2, "cooks")
)

expect_silent(
  plot(Model2, "hatplot")
)

expect_silent(
  plot(Model3, "qq")
)

expect_silent(
  plot(Model3, "marginal")
)

expect_silent(
  plot(Model3, "fitresid")
)

expect_error(
  plot(Model3, "bootHist")
)

expect_silent(
  plot(Model3, "rootogram")
)

expect_silent(
  plot(Model3, "dfpopContr")
)

expect_silent(
  plot(Model3, "dfpopBox")
)

expect_silent(
  plot(Model3, "scaleLoc")
)

expect_silent(
  plot(Model3, "cooks")
)

expect_silent(
  plot(Model3, "hatplot")
)

expect_silent(
  plot(Model4, "qq")
)

expect_silent(
  plot(Model4, "marginal")
)

expect_silent(
  plot(Model4, "fitresid")
)

expect_error(
  plot(Model4, "bootHist")
)

expect_silent(
  plot(Model4, "rootogram")
)

expect_silent(
  plot(Model4, "scaleLoc")
)

expect_error(
  plot(Model4, "cooks")
)

expect_silent(
  plot(Model4, "hatplot")
)

expect_silent(
  plot(Model5, "qq")
)

expect_silent(
  plot(Model5, "marginal")
)

expect_silent(
  plot(Model5, "fitresid")
)

expect_error(
  plot(Model5, "bootHist")
)

expect_silent(
  plot(Model5, "rootogram")
)

expect_silent(
  plot(Model5, "scaleLoc")
)

expect_error(
  plot(Model5, "cooks")
)

expect_silent(
  plot(Model5, "hatplot")
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
  tolerance = .5
)
