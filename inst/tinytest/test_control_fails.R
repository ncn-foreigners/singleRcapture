expect_error(
  controlPopVar(
    alpha = -.5
  )
)

expect_error(
  controlPopVar(
    alpha = c(.2, .5)
  )
)

expect_error(
  controlPopVar(
    sd = "abc"
  )
)

expect_error(
  controlPopVar(
    covType = "abc"
  )
)

expect_error(
  controlPopVar(
    bootType = "abc"
  )
)

expect_error(
  controlPopVar(
    B = -5
  )
)

expect_error(
  controlPopVar(
    traceBootstrapSize = -1
  )
)


expect_error(
  controlMethod(
    criterion = "abc"
  )
)

expect_error(
  controlMethod(
    criterion = c("abstol", "reltol")
  )
)

expect_error(
  controlMethod(
    epsilon = -1
  )
)

expect_error(
  controlMethod(
    maxiter = -1
  )
)

expect_error(
  controlMethod(
    verbose = c(2, 5)
  )
)

expect_silent(
  controlMethod(
    printEveryN = 2.5
  )
)

expect_error(
  controlMethod(
    etaStart = "abc"
  )
)

expect_error(
  controlMethod(
    coefStart = "abc"
  )
)

expect_error(
  controlMethod(
    silent = "abc"
  )
)

expect_error(
  controlMethod(
    checkDiagWeights = "abc"
  )
)

expect_error(
  controlMethod(
    saveIRLSlogs = "abc"
  )
)

expect_error(
  controlMethod(
    stepsize = -1
  )
)

expect_error(
  controlMethod(
    momentumFactor = -5
  )
)

expect_error(
  controlMethod(
    momentumActivation = -1
  )
)
