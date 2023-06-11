# cheking popsize estimation as in van der heijden et.al (2018)
set.seed(123)

expect_silent(
  opt <- estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "ztpoisson",
    method = "optim",
    popVar = "analytic",
    data = farmsubmission,
    controlMethod = controlMethod(
      maxiter = 5000,
      optimMethod = "Nelder-Mead",
      silent = TRUE
    )
  )
)

expect_silent(
  irls <- estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "ztpoisson",
    method = "IRLS",
    popVar = "analytic",
    data = farmsubmission,
    controlMethod = controlMethod(
      silent = TRUE
    )
  )
)

expect_equal(
  round(irls$populationSize$pointEstimate,
        digits = 0),
  18346
)

expect_equal(
  round(opt$populationSize$pointEstimate,
        digits = 0),
  18346
)
expect_equivalent(
  c(17932, 18760),
  as.numeric(round(irls$populationSize$confidenceInterval[1, ],
        digits = 0))
)
expect_equivalent(
  c(17932, 18760),
  as.numeric(round(opt$populationSize$confidenceInterval[1, ],
                   digits = 0))
)

expect_silent(
  ch <- estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "chao",
    method = "IRLS",
    popVar = "analytic",
    data = farmsubmission,
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_equal(
  round(ch$populationSize$pointEstimate,
        digits = 0),
  21657
)
expect_equivalent(
  c(20883, 22430),
  as.numeric(round(ch$populationSize$confidenceInterval[1, ],
                   digits = 0))
)

# on netherlandsimmigrant
# confint

expect_silent(
  zl <- estimatePopsize(
    formula = capture ~ . - reason - nation,
    data = netherlandsimmigrant,
    model = "zelterman",
    popVar = "analytic",
    method = "IRLS",
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_silent(
  zl2 <- estimatePopsize(
    formula = capture ~ 1,
    data = netherlandsimmigrant,
    model = "zelterman",
    popVar = "analytic",
    method = "IRLS",
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_silent(
  poi <- estimatePopsize(
    formula = capture ~ .,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "analytic",
    method = "IRLS",
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_silent(
  poi2 <- estimatePopsize(
    formula = capture ~ 1,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "analytic",
    method = "IRLS",
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_equivalent(
  as.numeric(round(zl$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8416, 12009)
)

expect_equivalent(
  as.numeric(round(zl2$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(8084, 10765)
)

expect_equivalent(
  as.numeric(round(poi$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(7185, 18198)
)

expect_equivalent(
  as.numeric(round(poi2$populationSize$confidenceInterval[1, ],
                   digits = 0)),
  c(6363, 7797)
)

# point
expect_equal(
  round(poi$populationSize$pointEstimate,
        digits = 0),
  12691
)

expect_equal(
  round(poi2$populationSize$pointEstimate,
        digits = 0),
  7080
)

expect_equal(
  round(zl$populationSize$pointEstimate,
        digits = 0),
  10213
)

expect_equal(
  round(zl2$populationSize$pointEstimate,
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
    controlMethod = controlMethod(epsilon = 1e-6, silent = TRUE),# testing silent
    controlPopVar = controlPopVar(B = 10, bootstrapFitcontrol = controlMethod(silent = TRUE, epsilon = .Machine$double.eps))
  )
)

expect_silent(
  estimatePopsize(
    formula = TOTAL_SUB ~ .,
    data = farmsubmission,
    model = "ztoigeom",
    popVar = "bootstrap",
    method = "IRLS",
    controlMethod = controlMethod(epsilon = 1e-6, silent = TRUE),
    controlPopVar = controlPopVar(
      B = 10, 
      bootstrapFitcontrol = controlMethod(
        silent = TRUE, 
        epsilon = .Machine$double.eps
      )
    )
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
    formula = capture ~ gender,
    data = netherlandsimmigrant,
    model = "ztpoisson",
    popVar = "bootstrap",
    method = "IRLS",
    controlPopVar = controlPopVar(
      B = 10, 
      bootstrapFitcontrol = controlMethod(
        silent = TRUE, 
        epsilon = .Machine$double.eps)
      )
  )
)

expect_error(
  estimatePopsize(formula = cbind(TOTAL_SUB, log_size) ~ log_distance, 
                  data = farmsubmission, model = ztpoisson)
)
