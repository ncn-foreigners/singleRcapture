# These tests are only supposed to be run on developer's machine and 
# package GitHub page not on CRAN (they take too long)

if (isTRUE(tolower(Sys.getenv("TEST_SINGLERCAPTURE_MULTICORE_DEVELOPER")) == "true")) {
  expect_error(
    estimatePopsize(
      formula = TOTAL_SUB ~ .,
      data = farmsubmission,
      model = "zotpoisson",
      method = "maxLik"
    )
  )
  
  expect_silent(
    xx <- estimatePopsize(
      formula = TOTAL_SUB ~ .,
      data = farmsubmission,
      model = "zotpoisson",
      controlMethod = controlMethod(
        epsilon = 1e-6, silent = TRUE
      )# testing silent
    )
  )
  
  expect_silent(
    summary(marginalFreq(xx), dropl5 = "group")
  )
  
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
        cores = 2L,
        confType = "normal"
      )
    )
  )
  
  expect_silent(
    predict(
      xx,
      type = "mean",
      se.fit = TRUE
    )
  )
  
  expect_warning(
    xx <- estimatePopsize(
      formula = TOTAL_SUB ~ log_size,
      data = farmsubmission,
      model = "zotnegbin",
      controlModel  = controlModel(alphaFormula = ~ C_TYPE), 
      controlMethod = controlMethod(
        silent = TRUE, stepsize = .7, 
        verbose = 4, maxiter = 20
      )
    )
  )
  
  expect_silent(
    predict(
      xx,
      type = "mean",
      se.fit = TRUE
    )
  )
  
  expect_silent(
    predict(
      xx,
      type = "response",
      se.fit = TRUE
    )
  )
  
  expect_silent(
    predict(
      xx,
      type = "link",
      se.fit = TRUE
    )
  )
  
  set.seed(123)
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
        bootType = "nonparametric",
        confType = "basic"
      )
    )
  )
  
  
  expect_silent(
    dfbeta(xx)
  )
  
  expect_silent(
    dfbeta(xx, cores = 2L)
  )
  
  expect_silent(
    predict(
      xx,
      type = "mean",
      se.fit = TRUE
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
  
  expect_warning(
    xx <- estimatePopsize(
      formula = TOTAL_SUB ~ C_TYPE, 
      data = df, 
      model = ztpoisson,
      popVar = "bootstrap",
      weights = df$ww,
      controlMethod = controlMethod(
        verbose = 5, 
        saveIRLSlogs = TRUE,
        criterion = "reltol"
      ),
      controlModel = controlModel(weightsAsCounts = TRUE),
      controlPopVar = controlPopVar(
        B = 70,
        cores = 2L,
        bootType = "semiparametric"
      )
    )
  )
  
  expect_true(
    !is.null(
      xx$fittingLog
    )
  )
  
  expect_silent(
    dfbeta(xx)
  )
  
  expect_silent(
    dfbeta(xx, cores = 2L)
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