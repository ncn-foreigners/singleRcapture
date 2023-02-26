# cheking popsize estimation as in van der heijden et.al (2018)
set.seed(123)
expect_equal(
  round(estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "ztpoisson",
    method = "IRLS",
    popVar = "analytic",
    data = farmsubmission
  )$populationSize$pointEstimate,
        digits = 0),
  18346
)
expect_equal(
  round(estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "ztpoisson",
    method = "optim",
    popVar = "analytic",
    data = farmsubmission
  )$populationSize$pointEstimate,
        digits = 0),
  18346
)
expect_equivalent(
  c(17932, 18760),
  as.numeric(round(estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "ztpoisson",
    method = "IRLS",
    popVar = "analytic",
    data = farmsubmission
  )$populationSize$confidenceInterval[1, ],
        digits = 0))
)
expect_equivalent(
  c(17932, 18760),
  as.numeric(round(estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "ztpoisson",
    method = "optim",
    popVar = "analytic",
    data = farmsubmission
  )$populationSize$confidenceInterval[1, ],
                   digits = 0))
)
expect_equal(
  round(estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "chao",
    method = "IRLS",
    popVar = "analytic",
    data = farmsubmission,
    controlMethod = controlMethod(silent = TRUE)
  )$populationSize$pointEstimate,
        digits = 0),
  21657
)
expect_equivalent(
  c(20883, 22430),
  as.numeric(round(estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "chao",
    method = "IRLS",
    popVar = "analytic",
    data = farmsubmission,
    controlMethod = controlMethod(silent = TRUE)
  )$populationSize$confidenceInterval[1, ],
                   digits = 0))
)
expect_equal(
  round(estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "chao",
    method = "optim",
    popVar = "analytic",
    data = farmsubmission
  )$populationSize$pointEstimate,
        digits = 0),
  21657
)
expect_equivalent(
  c(20883, 22430),
  as.numeric(round(estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "chao",
    method = "optim",
    popVar = "analytic",
    data = farmsubmission
  )$populationSize$confidenceInterval[1, ],
                   digits = 0))
)
# on netherlands imigrant
# confint
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ .,
    data = netherlandsimmigrant,
    model = "zelterman",
    popVar = "analytic",
    method = "IRLS"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(9983, 22394)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ .,
    data = netherlandsimmigrant,
    model = "zelterman",
    popVar = "analytic",
    method = "optim"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(9983, 22394)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ . - reason,
    data = netherlandsimmigrant,
    model = "zelterman",
    popVar = "analytic",
    method = "IRLS"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(9973, 22285)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ . - reason,
    data = netherlandsimmigrant,
    model = "zelterman",
    popVar = "analytic",
    method = "optim"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(9973, 22285)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ . - reason - nation,
    data = netherlandsimmigrant,
    model = "zelterman",
    popVar = "analytic",
    method = "IRLS",
    controlMethod = controlMethod(silent = TRUE)
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8416, 12009)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ . - reason - nation,
    data = netherlandsimmigrant,
    model = "zelterman",
    popVar = "analytic",
    method = "optim"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8416, 12009)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ gender,
    data = netherlandsimmigrant,
    model = "zelterman",
    popVar = "analytic",
    method = "IRLS"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8327, 11614)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ gender,
    data = netherlandsimmigrant,
    model = "zelterman",
    popVar = "analytic",
    method = "optim"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8327, 11614)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ 1,
    data = netherlandsimmigrant,
    model = "zelterman",
    popVar = "analytic",
    method = "IRLS"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8084, 10765)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ 1,
    data = netherlandsimmigrant,
    model = "zelterman",
    popVar = "analytic",
    method = "optim"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8084, 10765)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ .,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "analytic",
    method = "IRLS"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(7185, 18198)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ .,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "analytic",
    method = "optim"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(7185, 18198)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ . - reason,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "analytic",
    method = "IRLS",
    controlMethod = controlMethod(silent = TRUE)
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
c(7186, 18194)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ . - reason,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "analytic",
    method = "optim"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(7186, 18194)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ . - reason - nation,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "analytic",
    method = "IRLS",
    controlMethod = controlMethod(silent = TRUE)
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6637, 8978)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ gender,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "analytic",
    method = "optim"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6504, 8134)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ gender,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "analytic",
    method = "IRLS"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6504, 8134)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ gender,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "analytic",
    method = "optim"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6504, 8134)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ 1,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "analytic",
    method = "IRLS"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6363, 7797)
)
expect_equivalent(
  as.numeric(round(estimatePopsize(
    formula = capture ~ 1,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "analytic",
    method = "optim"
  )$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6363, 7797)
)
# point
expect_equal(
  round(estimatePopsize(
    formula = capture ~ .,
    model = "ztpoisson",
    method = "IRLS",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  12691
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ .,
    model = "ztpoisson",
    method = "optim",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  12691
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ . - reason,
    model = "ztpoisson",
    method = "IRLS",
    popVar = "analytic",
    data = netherlandsimmigrant,
    controlMethod = controlMethod(silent = TRUE)
  )$populationSize$pointEstimate,
        digits = 0),
  12690
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ . - reason,
    model = "ztpoisson",
    method = "optim",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  12690
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ . - reason - nation,
    model = "ztpoisson",
    method = "IRLS",
    popVar = "analytic",
    data = netherlandsimmigrant,
    controlMethod = controlMethod(silent = TRUE)
  )$populationSize$pointEstimate,
        digits = 0),
  7807
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ .- reason - nation,
    model = "ztpoisson",
    method = "optim",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  7807
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ gender,
    model = "ztpoisson",
    method = "IRLS",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  7319
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ gender,
    model = "ztpoisson",
    method = "optim",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  7319
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ 1,
    model = "ztpoisson",
    method = "IRLS",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  7080
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ 1,
    model = "ztpoisson",
    method = "optim",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  7080
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ .,
    model = "zelterman",
    method = "IRLS",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  16188
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ .,
    model = "zelterman",
    method = "optim",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  16188
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ . - reason,
    model = "zelterman",
    method = "IRLS",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  16129
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ . - reason,
    model = "zelterman",
    method = "optim",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  16129
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ . - reason - nation,
    model = "zelterman",
    method = "IRLS",
    popVar = "analytic",
    data = netherlandsimmigrant,
    controlMethod = controlMethod(silent = TRUE)
  )$populationSize$pointEstimate,
        digits = 0),
  10213
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ .- reason - nation,
    model = "zelterman",
    method = "optim",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  10213
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ gender,
    model = "zelterman",
    method = "IRLS",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  9970
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ gender,
    model = "zelterman",
    method = "optim",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  9970
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ 1,
    model = "zelterman",
    method = "IRLS",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  9425
)
expect_equal(
  round(estimatePopsize(
    formula = capture ~ 1,
    model = "zelterman",
    method = "optim",
    popVar = "analytic",
    data = netherlandsimmigrant
  )$populationSize$pointEstimate,
        digits = 0),
  9425
)
# other tests
expect_error(
  estimatePopsize(
    formula = Y ~ X^2, 
    data = data.frame(Y = rep(c(0, 1, 2, 3), 50),
                      X = rep(c(1, 2), 100)),
    model = "ztnegbin",
    method = "IRLS",
    popVar = "analytic"
  )
)
expect_silent(
  estimatePopsize(
    formula = TOTAL_SUB ~ .,
    data = farmsubmission,
    model = "zotpoisson",
    popVar = "bootstrap",
    method = "IRLS",
    controlMethod = controlMethod(epsilon = 1e-6),
    controlPopVar = controlPopVar(B = 10, bootstrapFitcontrol = controlMethod(silent = TRUE, epsilon = .Machine$double.eps))
  )
)
expect_silent(
  estimatePopsize(
    formula = TOTAL_SUB ~ .,
    data = farmsubmission,
    model = "ztnegbin",
    popVar = "bootstrap",
    method = "optim",
    controlMethod = controlMethod(epsilon = 1e-6),
    controlPopVar = controlPopVar(B = 10, bootstrapFitcontrol = controlMethod(silent = TRUE, epsilon = .Machine$double.eps))
  )
)
expect_silent(
  estimatePopsize(
    formula = TOTAL_SUB ~ .,
    data = farmsubmission,
    model = "chao",
    popVar = "bootstrap",
    method = "IRLS",
    controlMethod = controlMethod(silent = TRUE),
    controlPopVar = controlPopVar(B = 10, bootstrapFitcontrol = controlMethod(silent = TRUE, epsilon = .Machine$double.eps))
  )
)
expect_silent(
  estimatePopsize(
    formula = TOTAL_SUB ~ .,
    data = farmsubmission,
    model = "zelterman",
    popVar = "bootstrap",
    method = "IRLS",
    controlMethod = controlMethod(epsilon = 1e-6, silent = TRUE),
    controlPopVar = controlPopVar(B = 10, bootstrapFitcontrol = controlMethod(silent = TRUE, epsilon = .Machine$double.eps))
  )
)
expect_silent(
  estimatePopsize(
    formula = capture ~ gender,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "bootstrap",
    method = "IRLS",
    controlPopVar = controlPopVar(B = 10, bootstrapFitcontrol = controlMethod(silent = TRUE, epsilon = .Machine$double.eps))
  )
)

expect_error(
  estimatePopsize(formula = cbind(TOTAL_SUB, log_size) ~ log_distance, 
                   data = farmsubmission, model = ztpoisson)
)