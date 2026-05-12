# chao family tests

chao_family <- chao()
chao_eta <- matrix(chao_family$links[[1]](1.7))

expect_equal(
  as.numeric(chao_family$mu.eta(chao_eta, type = "trunc")),
  1.7 / (2 + 1.7),
  tolerance = 1e-10
)

expect_equal(
  as.numeric(chao_family$mu.eta(chao_eta, type = "nontrunc")),
  1.7,
  tolerance = 1e-10
)

expect_equal(
  as.numeric(chao_family$variance(chao_eta, type = "trunc")),
  (1.7 / (2 + 1.7)) * (2 / (2 + 1.7)),
  tolerance = 1e-10
)

expect_equal(
  as.numeric(chao_family$densityFunction(1, chao_eta, type = "trunc")),
  dpois(1, 1.7) / (1 - dpois(0, 1.7)),
  tolerance = 1e-10
)

expect_equal(
  as.numeric(chao_family$densityFunction(2, chao_eta, type = "nontrunc")),
  dpois(2, 1.7),
  tolerance = 1e-10
)

chao_dev_y2 <- as.numeric(chao_family$devResids(
  y = 2, eta = matrix(chao_family$links[[1]](100)), wt = 1
))
chao_dev_y1 <- as.numeric(chao_family$devResids(
  y = 1, eta = matrix(chao_family$links[[1]](0.01)), wt = 1
))
expect_true(is.finite(chao_dev_y2))
expect_true(is.finite(chao_dev_y1))
expect_true(chao_dev_y2 > 0)
expect_true(chao_dev_y1 < 0)
expect_true(abs(chao_dev_y2) < 1)
expect_true(abs(chao_dev_y1) < 1)

expect_equal(
  as.numeric(chao_family$devResids(y = 0, eta = chao_eta, wt = 1)),
  0,
  tolerance = 1e-10
)
expect_equal(
  as.numeric(chao_family$devResids(y = 3, eta = chao_eta, wt = 1)),
  0,
  tolerance = 1e-10
)

expect_silent(
  fit_chao <- estimatePopsize(
    formula = TOTAL_SUB ~ .,
    model = chao(),
    data = farmsubmission,
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_silent(
  predict(
    fit_chao,
    type = "mean"
  )
)

expect_equivalent(
  fit_chao$populationSize$pointEstimate,
  21657,
  tolerance = .005
)

expect_equivalent(
  fit_chao$populationSize$confidenceInterval[1, ],
  list(20883, 22430),
  tolerance = .005
)

expect_silent(
  population_chao <- popSizeEst(fit_chao)
)

expect_true(is.finite(population_chao$variance))
expect_true(population_chao$variance > 0)
