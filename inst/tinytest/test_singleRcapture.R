# cheking popsize estimation as in van der heijden et.al (2018)

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
  estimatePopsize(
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

# These tests are only supposed to be run on developer's machine and 
# package GitHub page not on CRAN (they take too long)

if (isTRUE(tolower(Sys.getenv("TEST_SINGLERCAPTURE_MULTICORE_DEVELOPER")) == "true")) {
  set.seed(123)
  expect_silent(
    estimatePopsize(
      formula = TOTAL_SUB ~ .,
      data = farmsubmission,
      model = "zotpoisson",
      popVar = "bootstrap",
      controlMethod = controlMethod(epsilon = 1e-6, silent = TRUE),# testing silent
      controlPopVar = controlPopVar(
        B = 140, 
        bootstrapFitcontrol = controlMethod(
          silent = TRUE, 
          epsilon = .Machine$double.eps
        ),
        cores = 2L
      )
    )
  )
  
  expect_silent(
    estimatePopsize(
      formula = capture ~ gender,
      data = netherlandsimmigrant,
      model = "zotgeom",
      method = "optim",
      popVar = "bootstrap",
      controlMethod = controlMethod(epsilon = 1e-6, silent = TRUE),# testing silent
      controlPopVar = controlPopVar(
        B = 70,
        cores = 2L, bootstrapFitcontrol = controlMethod(),
        bootType = "nonparametric"
      )
    )
  )
  
  expect_silent(
    estimatePopsize(
      formula = TOTAL_SUB ~ .,
      data = farmsubmission,
      model = "oiztgeom",
      popVar = "bootstrap",
      controlModel = controlModel(omegaFormula = ~ .),
      controlMethod = controlMethod(silent = TRUE),# testing silent
      controlPopVar = controlPopVar(
        B = 35,
        cores = 2L, bootstrapFitcontrol = controlMethod(),
        bootType = "semiparametric"
      )
    )
  )
  
  df <- farmsubmission[, c(1,4)]
  df$ww <- 0
  ### this is dplyr::count but slower and without dependencies
  df <- aggregate(ww ~ ., df, FUN = length)
  
  expect_silent(
    estimatePopsize(
      formula = TOTAL_SUB ~ C_TYPE, 
      data = df, 
      model = ztpoisson,
      popVar = "bootstrap",
      weights = df$ww,
      controlMethod = controlMethod(silent = TRUE),
      controlModel = controlModel(weightsAsCounts = TRUE),
      controlPopVar = controlPopVar(
        B = 70,
        cores = 2L,
        bootType = "semiparametric"
      )
    )
  )
  
  expect_silent(
    estimatePopsize(
      formula = TOTAL_SUB ~ C_TYPE, 
      data = df, 
      model = ztoipoisson,
      popVar = "bootstrap",
      weights = df$ww,
      controlMethod = controlMethod(silent = TRUE),
      controlModel = controlModel(weightsAsCounts = TRUE),
      controlPopVar = controlPopVar(
        B = 70,
        cores = 2L,
        bootType = "parametric"
      )
    )
  )
  
  expect_silent(
    estimatePopsize(
      formula = TOTAL_SUB ~ C_TYPE, 
      data = df, 
      model = ztoigeom,
      popVar = "bootstrap",
      weights = df$ww,
      controlMethod = controlMethod(silent = TRUE),
      controlModel = controlModel(weightsAsCounts = TRUE),
      controlPopVar = controlPopVar(
        B = 70,
        cores = 2L,
        bootType = "nonparametric"
      )
    )
  )
  
  expect_silent(
    estimatePopsize(
      formula = TOTAL_SUB ~ C_TYPE, 
      data = df, 
      model = chao,
      popVar = "bootstrap",
      weights = df$ww,
      controlMethod = controlMethod(silent = TRUE),
      controlModel = controlModel(weightsAsCounts = TRUE),
      controlPopVar = controlPopVar(
        B = 70,
        bootType = "nonparametric"
      )
    )
  )
  
  expect_silent(
    estimatePopsize(
      formula = TOTAL_SUB ~ C_TYPE, 
      data = df, 
      model = zelterman,
      popVar = "bootstrap",
      weights = df$ww,
      controlMethod = controlMethod(silent = TRUE),
      controlModel = controlModel(weightsAsCounts = TRUE),
      controlPopVar = controlPopVar(
        B = 70,
        bootType = "semiparametric"
      )
    )
  )
  
  expect_silent(
    estimatePopsize(
      formula = TOTAL_SUB ~ C_TYPE, 
      data = df, 
      model = ztgeom,
      popVar = "bootstrap",
      weights = df$ww,
      controlMethod = controlMethod(silent = TRUE),
      controlModel = controlModel(weightsAsCounts = TRUE),
      controlPopVar = controlPopVar(
        B = 70,
        bootType = "parametric"
      )
    )
  )
}