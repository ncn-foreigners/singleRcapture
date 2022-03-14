# cheking popsize estimation as in van der heijden et.al (2018)
load("~/Desktop/singleRcapture/inst/tinytest/carcassubmission.rda")
load("~/Desktop/singleRcapture/inst/tinytest/netherlandsimmigrant.rda")
chaocovtotal <- carcassubmisssion
dutchc <- netherlandsimmigrant
expect_equal(
  round(estimate_popsize(formula = TOTAL_SUB ~ .,
                         model = "ztpoisson",
                         method = "robust",
                         pop.var = "analytic",
                         data = chaocovtotal)$populationSize$pointEstimate,
        digits = 0),
  18346
)
expect_equivalent(
  c(17932, 18760),
  round(estimate_popsize(formula = TOTAL_SUB ~ .,
                   model = "ztpoisson",
                   method = "mle",
                   pop.var = "analytic",
                   data = chaocovtotal)$populationSize$confidenceInterval,
        digits = 0)
)
expect_equivalent(
  round(estimate_popsize(formula = capture ~ .,
                         data = dutchc,
                         model = "zelterman",
                         pop.var = "analytic",
                         method = "mle")$populationSize$confidenceInterval,
        digits = 0),
  c(9983, 22394)
)
