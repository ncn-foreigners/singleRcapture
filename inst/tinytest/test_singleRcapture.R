# cheking popsize estimation as in van der heijden et.al (2018)
expect_equal(
  round(estimate_popsize(formula = TOTAL_SUB ~ .,
                         model = "ztpoisson",
                         method = "robust",
                         pop.var = "analytic",
                         data = farmsubmission)$populationSize$pointEstimate,
        digits = 0),
  18346
)
expect_equal(
  round(estimate_popsize(formula = TOTAL_SUB ~ .,
                         model = "ztpoisson",
                         method = "mle",
                         pop.var = "analytic",
                         data = farmsubmission)$populationSize$pointEstimate,
        digits = 0),
  18346
)
expect_equivalent(
  c(17932, 18760),
  as.numeric(round(estimate_popsize(formula = TOTAL_SUB ~ .,
                   model = "ztpoisson",
                   method = "robust",
                   pop.var = "analytic",
                   data = farmsubmission)$populationSize$confidenceInterval[1, ],
        digits = 0))
)
expect_equivalent(
  c(17932, 18760),
  as.numeric(round(estimate_popsize(formula = TOTAL_SUB ~ .,
                                    model = "ztpoisson",
                                    method = "mle",
                                    pop.var = "analytic",
                                    data = farmsubmission)$populationSize$confidenceInterval[1, ],
                   digits = 0))
)
expect_equivalent(
  as.numeric(round(estimate_popsize(formula = capture ~ .,
                         data = netherlandsimmigrant,
                         model = "zelterman",
                         pop.var = "analytic",
                         method = "robust")$populationSize$confidenceInterval[1, ],
        digits = 0)),
  c(9983, 22394)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(formula = capture ~ .,
                                    data = netherlandsimmigrant,
                                    model = "zelterman",
                                    pop.var = "analytic",
                                    method = "mle")$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(9983, 22394)
)
expect_error(
  estimate_popsize(formula = Y ~ X^2, 
                   data = data.frame(Y = rep(c(0, 1, 2, 3), 50),
                                     X = rep(c(1, 2), 100)),
                   model = "ztnegbin",
                   method = "robust",
                   pop.var = "analytic")
)
expect_silent(
  estimate_popsize(formula = TOTAL_SUB ~ .,
                   data = farmsubmission,
                   model = "zotpoisson",
                   pop.var = "bootstrap",
                   method = "robust",
                   control.pop.var = control.pop.var(strapNumber = 10))
)
expect_silent(
  estimate_popsize(formula = TOTAL_SUB ~ .,
                   data = farmsubmission,
                   model = "ztnegbin",
                   pop.var = "bootstrap",
                   method = "mle",
                   control.pop.var = control.pop.var(strapNumber = 10))
)
expect_silent(
  estimate_popsize(formula = TOTAL_SUB ~ .,
                   data = farmsubmission,
                   model = "chao",
                   pop.var = "bootstrap",
                   method = "robust",
                   control.pop.var = control.pop.var(strapNumber = 10))
)
expect_silent(
  estimate_popsize(formula = TOTAL_SUB ~ .,
                   data = farmsubmission,
                   model = "zelterman",
                   pop.var = "bootstrap",
                   method = "robust",
                   control.pop.var = control.pop.var(strapNumber = 10))
)
expect_silent(
  estimate_popsize(formula = capture ~ .,
                   data = netherlandsimmigrant,
                   model = "ztpoisson",
                   pop.var = "bootstrap",
                   method = "robust",
                   control.pop.var = control.pop.var(strapNumber = 10))
)