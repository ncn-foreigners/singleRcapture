# These tests are only supposed to be run on developer's machine and 
# package GitHub page not on CRAN (they take too long)

if (isTRUE(tolower(Sys.getenv("TEST_SINGLERCAPTURE_MULTICORE_DEVELOPER")) == "true")) {
#if (TRUE) {
  set.seed(123)
  expect_silent(
    xx <- estimatePopsize(
      formula = TOTAL_SUB ~ .,
      data = farmsubmission,
      model = "zotpoisson",
      popVar = "bootstrap",
      controlMethod = controlMethod(epsilon = 1e-6, silent = TRUE),# testing silent
      controlPopVar = controlPopVar(
        B = 140, 
        bootType = "parametric",
        bootstrapFitcontrol = controlMethod(
          silent = TRUE, 
          epsilon = .Machine$double.eps
        ),
        cores = 2L
      )
    )
  )
  
  expect_silent(
    predict(
      xx,
      type = "mean"
    )
  )
  
  expect_silent(
    xx <- estimatePopsize(
      formula = capture ~ gender,
      data = netherlandsimmigrant,
      model = "zotgeom",
      method = "optim",
      popVar = "bootstrap",
      controlMethod = controlMethod(epsilon = 1e-6, silent = TRUE),# testing silent
      controlPopVar = controlPopVar(
        B = 70,
        cores = 2L, 
        bootstrapFitcontrol = controlMethod(),
        bootType = "nonparametric"
      )
    )
  )
  
  expect_silent(
    predict(
      xx,
      type = "mean"
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
        cores = 2L, 
        bootstrapFitcontrol = controlMethod(),
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