# cheking popsize estimation as in van der heijden et.al (2018)
set.seed(123)
expect_equal(
  round(estimate_popsize(
    formula = TOTAL_SUB ~ .,
    model = "ztpoisson",
    method = "robust",
    pop.var = "analytic",
    data = farmsubmission
  )$populationSize$pointEstimate,
        digits = 0),
  18346
)
expect_equal(
  round(estimate_popsize(
    formula = TOTAL_SUB ~ .,
    model = "ztpoisson",
    method = "mle",
    pop.var = "analytic",
    data = farmsubmission
  )$populationSize$pointEstimate,
        digits = 0),
  18346
)
expect_equivalent(
  c(17932, 18760),
  as.numeric(round(estimate_popsize(
    formula = TOTAL_SUB ~ .,
    model = "ztpoisson",
    method = "robust",
    pop.var = "analytic",
    data = farmsubmission
  )$populationSize$confidenceInterval[1, ],
        digits = 0))
)
expect_equivalent(
  c(17932, 18760),
  as.numeric(round(estimate_popsize(
    formula = TOTAL_SUB ~ .,
    model = "ztpoisson",
    method = "mle",
    pop.var = "analytic",
    data = farmsubmission
  )$populationSize$confidenceInterval[1, ],
                   digits = 0))
)
expect_equal(
  round(estimate_popsize(
    formula = TOTAL_SUB ~ .,
    model = "chao",
    method = "robust",
    pop.var = "analytic",
    data = farmsubmission,
    control.method = control.method(silent = TRUE)
  )$populationSize$pointEstimate,
        digits = 0),
  21657
)
expect_equivalent(
  c(20883, 22430),
  as.numeric(round(estimate_popsize(
    formula = TOTAL_SUB ~ .,
    model = "chao",
    method = "robust",
    pop.var = "analytic",
    data = farmsubmission,
    control.method = control.method(silent = TRUE)
  )$populationSize$confidenceInterval[1, ],
                   digits = 0))
)
expect_equal(
  round(estimate_popsize(
    formula = TOTAL_SUB ~ .,
    model = "chao",
    method = "mle",
    pop.var = "analytic",
    data = farmsubmission
  )$populationSize$pointEstimate,
        digits = 0),
  21657
)
expect_equivalent(
  c(20883, 22430),
  as.numeric(round(estimate_popsize(
    formula = TOTAL_SUB ~ .,
    model = "chao",
    method = "mle",
    pop.var = "analytic",
    data = farmsubmission
  )$populationSize$confidenceInterval[1, ],
                   digits = 0))
)
# on netherlands imigrant
# confint
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ .,
    data = netherlandsimmigrant,
    model = "zelterman",
    pop.var = "analytic",
    method = "robust"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(9983, 22394)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ .,
    data = netherlandsimmigrant,
    model = "zelterman",
    pop.var = "analytic",
    method = "mle"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(9983, 22394)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ . - reason,
    data = netherlandsimmigrant,
    model = "zelterman",
    pop.var = "analytic",
    method = "robust"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(9973, 22285)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ . - reason,
    data = netherlandsimmigrant,
    model = "zelterman",
    pop.var = "analytic",
    method = "mle"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(9973, 22285)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ . - reason - nation,
    data = netherlandsimmigrant,
    model = "zelterman",
    pop.var = "analytic",
    method = "robust",
    control.method = control.method(silent = TRUE)
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8416, 12009)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ . - reason - nation,
    data = netherlandsimmigrant,
    model = "zelterman",
    pop.var = "analytic",
    method = "mle"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8416, 12009)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ gender,
    data = netherlandsimmigrant,
    model = "zelterman",
    pop.var = "analytic",
    method = "robust"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8327, 11614)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ gender,
    data = netherlandsimmigrant,
    model = "zelterman",
    pop.var = "analytic",
    method = "mle"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8327, 11614)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ 1,
    data = netherlandsimmigrant,
    model = "zelterman",
    pop.var = "analytic",
    method = "robust"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8084, 10765)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ 1,
    data = netherlandsimmigrant,
    model = "zelterman",
    pop.var = "analytic",
    method = "mle"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8084, 10765)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ .,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    pop.var = "analytic",
    method = "robust"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(7185, 18198)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ .,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    pop.var = "analytic",
    method = "mle"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(7185, 18198)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ . - reason,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    pop.var = "analytic",
    method = "robust",
    control.method = control.method(silent = TRUE)
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
c(7186, 18194)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ . - reason,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    pop.var = "analytic",
    method = "mle"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(7186, 18194)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ . - reason - nation,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    pop.var = "analytic",
    method = "robust",
    control.method = control.method(silent = TRUE)
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6637, 8978)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ gender,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    pop.var = "analytic",
    method = "mle"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6504, 8134)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ gender,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    pop.var = "analytic",
    method = "robust"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6504, 8134)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ gender,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    pop.var = "analytic",
    method = "mle"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6504, 8134)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ 1,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    pop.var = "analytic",
    method = "robust"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6363, 7797)
)
expect_equivalent(
  as.numeric(round(estimate_popsize(
    formula = capture ~ 1,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    pop.var = "analytic",
    method = "mle"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6363, 7797)
)
# point
expect_equal(
  round(estimate_popsize(
    formula = capture ~ .,
    model = "ztpoisson",
    method = "robust",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  12691
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ .,
    model = "ztpoisson",
    method = "mle",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  12691
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ . - reason,
    model = "ztpoisson",
    method = "robust",
    pop.var = "analytic",
    data = netherlandsimmigrant,
    control.method = control.method(silent = TRUE)
  )$populationSize$pointEstimate,
        digits = 0),
  12690
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ . - reason,
    model = "ztpoisson",
    method = "mle",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  12690
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ . - reason - nation,
    model = "ztpoisson",
    method = "robust",
    pop.var = "analytic",
    data = netherlandsimmigrant,
    control.method = control.method(silent = TRUE)
  )$populationSize$pointEstimate,
        digits = 0),
  7807
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ .- reason - nation,
    model = "ztpoisson",
    method = "mle",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  7807
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ gender,
    model = "ztpoisson",
    method = "robust",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  7319
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ gender,
    model = "ztpoisson",
    method = "mle",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  7319
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ 1,
    model = "ztpoisson",
    method = "robust",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  7080
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ 1,
    model = "ztpoisson",
    method = "mle",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  7080
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ .,
    model = "zelterman",
    method = "robust",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  16188
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ .,
    model = "zelterman",
    method = "mle",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  16188
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ . - reason,
    model = "zelterman",
    method = "robust",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  16129
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ . - reason,
    model = "zelterman",
    method = "mle",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  16129
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ . - reason - nation,
    model = "zelterman",
    method = "robust",
    pop.var = "analytic",
    data = netherlandsimmigrant,
    control.method = control.method(silent = TRUE)
  )$populationSize$pointEstimate,
        digits = 0),
  10213
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ .- reason - nation,
    model = "zelterman",
    method = "mle",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  10213
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ gender,
    model = "zelterman",
    method = "robust",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  9970
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ gender,
    model = "zelterman",
    method = "mle",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  9970
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ 1,
    model = "zelterman",
    method = "robust",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  9425
)
expect_equal(
  round(estimate_popsize(
    formula = capture ~ 1,
    model = "zelterman",
    method = "mle",
    pop.var = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  9425
)
# other tests
expect_error(
  estimate_popsize(
    formula = Y ~ X^2, 
    data = data.frame(Y = rep(c(0, 1, 2, 3), 50),
                      X = rep(c(1, 2), 100)),
    model = "ztnegbin",
    method = "robust",
    pop.var = "analytic"
  )
)
expect_silent(
  estimate_popsize(
    formula = TOTAL_SUB ~ .,
    data = farmsubmission,
    model = "zotpoisson",
    pop.var = "bootstrap",
    method = "robust",
    control.method = control.method(epsilon = 1e-6),
    control.pop.var = control.pop.var(B = 10, bootstrapFitcontrol = control.method(silent = TRUE, epsilon = .Machine$double.eps))
  )
)
expect_silent(
  estimate_popsize(
    formula = TOTAL_SUB ~ .,
    data = farmsubmission,
    model = "ztnegbin",
    pop.var = "bootstrap",
    method = "mle",
    control.method = control.method(epsilon = 1e-6),
    control.pop.var = control.pop.var(B = 10, bootstrapFitcontrol = control.method(silent = TRUE, epsilon = .Machine$double.eps))
  )
)
expect_silent(
  estimate_popsize(
    formula = TOTAL_SUB ~ .,
    data = farmsubmission,
    model = "chao",
    pop.var = "bootstrap",
    method = "robust",
    control.method = control.method(silent = TRUE),
    control.pop.var = control.pop.var(B = 10, bootstrapFitcontrol = control.method(silent = TRUE, epsilon = .Machine$double.eps))
  )
)
expect_silent(
  estimate_popsize(
    formula = TOTAL_SUB ~ .,
    data = farmsubmission,
    model = "zelterman",
    pop.var = "bootstrap",
    method = "robust",
    control.method = control.method(epsilon = 1e-6, silent = TRUE),
    control.pop.var = control.pop.var(B = 10, bootstrapFitcontrol = control.method(silent = TRUE, epsilon = .Machine$double.eps))
  )
)
expect_silent(
  estimate_popsize(
    formula = capture ~ .,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    pop.var = "bootstrap",
    method = "robust",
    control.pop.var = control.pop.var(B = 10, bootstrapFitcontrol = control.method(silent = TRUE, epsilon = .Machine$double.eps))
  )
)