# cheking popsize estimation as in articles on which the package is based

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
  dd <- estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "zotnegbin",
    method = "optim",
    data = farmsubmission,
    controlMethod = controlMethod(
      maxiter = 5000,
      optimMethod = "Nelder-Mead",
      silent = TRUE
    )
  )
)

expect_silent(
  predict(
    dd,
    type = "mean"
  )
)

expect_silent(
  dd <- estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "zotgeom",
    method = "optim",
    data = farmsubmission
  )
)

expect_silent(
  predict(
    dd,
    type = "mean"
  )
)

expect_silent(
  irls <- estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = "ztpoisson",
    data = farmsubmission,
    controlMethod = controlMethod(
      silent = TRUE
    )
  )
)

expect_equivalent(
  irls$populationSize$pointEstimate,
  18346, tolerance = .005
)

expect_equivalent(
  opt$populationSize$pointEstimate,
  18346, tolerance = .005
)

expect_equivalent(
  irls$populationSize$confidenceInterval[1, ],
  list(17932, 18760),
  tolerance = .005
)

expect_equivalent(
  opt$populationSize$confidenceInterval[1, ],
  list(17932, 18760),
  tolerance = .005
)

expect_silent(
  ch <- estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = chao(),
    data = farmsubmission,
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_silent(
  predict(
    ch,
    type = "mean"
  )
)

expect_equivalent(
  ch$populationSize$pointEstimate,
  21657,
  tolerance = .005
)

expect_equivalent(
  ch$populationSize$confidenceInterval[1, ],
  list(20883, 22430),
  tolerance = .005
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
  predict(
    zl,
    type = "mean"
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
  zl$populationSize$confidenceInterval[1, ],
  list(8416, 12009),
  tolerance = .005
)

expect_equivalent(
  zl2$populationSize$confidenceInterval[1, ],
  list(8084, 10765),
  tolerance = .005
)

expect_equivalent(
  poi$populationSize$confidenceInterval[1, ],
  list(7185, 18198),
  tolerance = .005
)

expect_equivalent(
  poi2$populationSize$confidenceInterval[1, ],
  list(6363, 7797),
  tolerance = .005
)

# point
expect_equivalent(
  poi$populationSize$pointEstimate,
  12691,
  tolerance = .005
)

expect_equivalent(
  poi2$populationSize$pointEstimate,
  7080,
  tolerance = .005
)

expect_equivalent(
  zl$populationSize$pointEstimate,
  10213,
  tolerance = .005
)

expect_equivalent(
  zl2$populationSize$pointEstimate,
  9425,
  tolerance = .005
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

set.seed(123)
expect_silent(
  estimatePopsize(
    formula = TOTAL_SUB ~ .,
    data = farmsubmission,
    model = "zotpoisson",
    popVar = "bootstrap",
    controlMethod = controlMethod(epsilon = 1e-6, silent = TRUE),# testing silent
    controlPopVar = controlPopVar(B = 4, bootstrapFitcontrol = controlMethod(silent = TRUE, epsilon = .Machine$double.eps))
  )
)

expect_silent(
  estimatePopsize(
    formula = TOTAL_SUB ~ .,
    data = farmsubmission,
    model = "ztoigeom", 
    method = "optim",
    popVar = "bootstrap",
    controlMethod = controlMethod(epsilon = 1e-6, silent = TRUE),
    controlPopVar = controlPopVar(
      B = 6,
      bootType = "semiparametric"
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
    controlPopVar = controlPopVar(
      B = 6,
      bootType = "nonparametric"
    )
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
      B = 4, 
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

