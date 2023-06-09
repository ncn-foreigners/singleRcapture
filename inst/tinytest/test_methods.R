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

expect_silent(
  Model3 <- estimatePopsize(
    formula = capture ~ . - age, 
    data = netherlandsimmigrant, 
    model = chao, 
    method = "IRLS",
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_silent(
  Model4 <- estimatePopsize(
    formula = TOTAL_SUB ~ ., 
    data = farmsubmission, 
    model = ztoigeom, 
    method = "IRLS",
    controlPopVar = controlPopVar(covType = "Fisher"),
    controlModel = controlModel(omegaFormula = ~ log_distance),
    controlMethod = controlMethod(silent = TRUE)
  )
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

expect_silent(
  hatvalues(Model)
)

expect_silent(
  hatvalues(Model2)
)

expect_silent(
  hatvalues(Model3)
)

expect_silent(
  hatvalues(Model4)
)

expect_silent(
  dfp <- dfpopsize(Model, dfbeta = dfb)
)

expect_equal(
  max(abs(dfp)),
  4236.412,
  tolerance = .1
)

expect_equal(
  abs(mean(abs(dfp))),
  19.19,
  tolerance = .1
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

expect_true(
  abs(mean(abs(dfpopsize(Model2, dfbeta = dfb2))) - 16.59) < .1
)

expect_true(
  abs(max(abs(dfpopsize(Model3, dfbeta = dfb3))) - 3681.764) < .2
)

expect_true(
  abs(mean(abs(dfpopsize(Model3, dfbeta = dfb3))) - 17.18) < .1
)


# Extractors

expect_true(
  max(abs(
    c(AIC(Model), AIC(Model1), AIC(Model2), AIC(Model3), AIC(Model4)) -
      c(1712.901, 1805.904, 1131.723, 1133.006, 34580.82)
  )) < 1
)

expect_true(
  max(abs(
    c(BIC(Model), BIC(Model1), BIC(Model2), BIC(Model3), BIC(Model4)) -
      c(1757.213, 1811.443, 1170.3, 1177.094, 34625.2)
  )) < 1
)

expect_silent(
  c(extractAIC(Model), extractAIC(Model1), extractAIC(Model2), 
    extractAIC(Model3), extractAIC(Model4))
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

expect_false(
  all(vcov(Model, type = "observedInform") == vcov(Model, type = "Fisher"))
)

expect_false(
  all(vcov(Model1, type = "observedInform") == vcov(Model1, type = "Fisher"))
)

expect_false(
  all(vcov(Model4, type = "observedInform") == vcov(Model4, type = "Fisher"))
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
  max(bread(Model2, type = "observedInform") - bread(Model2, type = "Fisher")) < 1e-4
)

expect_true(
  max(bread(Model3, type = "observedInform") - bread(Model3, type = "Fisher")) < 1e-4
)

expect_false(
  all(bread(Model4, type = "observedInform") == bread(Model4, type = "Fisher"))
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

expect_identical(
  family(Model),
  Model$model
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
  all(dim(model.matrix(Model4, "vlm")) == c(2 * 12036, 6))
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

expect_equal(
  c(popSizeEst(Model)$pointEstimate, popSizeEst(Model1)$pointEstimate, 
    popSizeEst(Model2)$pointEstimate, popSizeEst(Model3)$pointEstimate,
    popSizeEst(Model4)$pointEstimate),
  c(12690, 7080, 15816.14, 15713.14, 29478.72),
  tol = .05
)

expect_equal(
  c(popSizeEst(Model)$variance, popSizeEst(Model1)$variance, 
    popSizeEst(Model2)$variance, popSizeEst(Model3)$variance,
    popSizeEst(Model4)$variance),
  c(7885812, 133774.1, 9093464, 9096077, 431426.9),
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
  summary(marginalFreq(Model3), df = 1, dropl5 = "group")
)

expect_true(
  all(summary(marginalFreq(Model3), df = 1, dropl5 = "group")$Test$`P(>X^2)` < .001)
)
