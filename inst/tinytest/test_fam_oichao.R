# oichao family tests

oichao_theory_data <- data.frame(
  y = c(rep(2, 95), rep(3, 32))
)

expect_silent(
  fit_oichao_theory <- estimatePopsize(
    formula = y ~ 1,
    model = "oichao",
    data = oichao_theory_data,
    method = "IRLS",
    popVar = "analytic",
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_equivalent(
  3 * exp(unname(coef(fit_oichao_theory))[1]),
  3 * 32 / 95,
  tolerance = .005
)

expect_equivalent(
  fit_oichao_theory$populationSize$pointEstimate,
  127 + (2 / 9) * 95 ^ 3 / 32 ^ 2,
  tolerance = .005
)

set.seed(2024)
data_oichao_covariate <- data.frame(
  x = rnorm(500)
)
data_oichao_covariate$lambda <- exp(0.2 + 0.35 * data_oichao_covariate$x)
data_oichao_covariate$y <- rpois(
  nrow(data_oichao_covariate),
  data_oichao_covariate$lambda
)
data_oichao_covariate <- subset(data_oichao_covariate, y > 0)

expect_silent(
  fit_oichao_irls <- estimatePopsize(
    formula = y ~ x,
    model = "oichao",
    data = data_oichao_covariate,
    method = "IRLS",
    popVar = "analytic",
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_silent(
  fit_oichao_optim <- estimatePopsize(
    formula = y ~ x,
    model = oichao(),
    data = data_oichao_covariate,
    method = "optim",
    popVar = "analytic",
    controlMethod = controlMethod(
      silent = TRUE,
      optimMethod = "BFGS",
      maxiter = 500
    )
  )
)

expect_true(
  max(abs(coef(fit_oichao_irls) - coef(fit_oichao_optim))) < 1e-3
)

expect_true(
  abs(
    fit_oichao_irls$populationSize$pointEstimate -
      fit_oichao_optim$populationSize$pointEstimate
  ) < 0.25
)

expect_silent(
  predict(
    fit_oichao_irls,
    type = "response",
    se.fit = TRUE
  )
)

expect_silent(summary(fit_oichao_irls))

expect_silent(popSizeEst(fit_oichao_irls))

set.seed(321)
observed_oichao <- sample(c(1, 2, 3, 4), 20, replace = TRUE)
design_oichao <- cbind(1, rnorm(20))
coefficients_oichao <- c(-0.3, 0.2)
family_oichao <- oichao()
objective_oichao <- family_oichao$makeMinusLogLike(
  y = observed_oichao,
  X = design_oichao,
  weight = rep(1, 20),
  deriv = 0
)
gradient_oichao <- family_oichao$makeMinusLogLike(
  y = observed_oichao,
  X = design_oichao,
  weight = rep(1, 20),
  deriv = 1
)
hessian_oichao <- family_oichao$makeMinusLogLike(
  y = observed_oichao,
  X = design_oichao,
  weight = rep(1, 20),
  deriv = 2
)
eps_fd_oichao <- 1e-6
finite_difference_gradient_oichao <- sapply(1:2, function(j) {
  coefficients_plus <- coefficients_oichao
  coefficients_minus <- coefficients_oichao
  coefficients_plus[j] <- coefficients_plus[j] + eps_fd_oichao
  coefficients_minus[j] <- coefficients_minus[j] - eps_fd_oichao
  (objective_oichao(coefficients_plus) - objective_oichao(coefficients_minus)) /
    (2 * eps_fd_oichao)
})
finite_difference_hessian_oichao <- matrix(0, 2, 2)
for (i in 1:2) {
  for (j in 1:2) {
    ei <- rep(0, 2)
    ej <- rep(0, 2)
    ei[i] <- eps_fd_oichao
    ej[j] <- eps_fd_oichao
    finite_difference_hessian_oichao[i, j] <- (
      objective_oichao(coefficients_oichao + ei + ej) -
        objective_oichao(coefficients_oichao + ei - ej) -
        objective_oichao(coefficients_oichao - ei + ej) +
        objective_oichao(coefficients_oichao - ei - ej)
    ) / (4 * eps_fd_oichao ^ 2)
  }
}

expect_true(
  max(abs(
    as.numeric(-gradient_oichao(coefficients_oichao)) -
      finite_difference_gradient_oichao
  )) < 1e-5
)

expect_true(
  max(abs(
    -hessian_oichao(coefficients_oichao) -
      finite_difference_hessian_oichao
  )) < 1e-3
)

expect_silent(
  fit_oichao <- estimatePopsize(
    formula = y ~ 1,
    data = data.frame(y = c(rep(2, 20), rep(3, 10))),
    model = "oichao",
    controlMethod = controlMethod(silent = TRUE)
  )
)

expect_true(
  all(as.matrix(simulate(fit_oichao, nsim = 4, seed = 1)) %in% 2:3)
)

eta_oichao <- matrix(family_oichao$links[[1]](1.7))

expect_equal(
  as.numeric(family_oichao$mu.eta(eta_oichao, type = "trunc")),
  1.7 / (3 + 1.7),
  tolerance = 1e-10
)

expect_equal(
  as.numeric(family_oichao$mu.eta(eta_oichao, type = "nontrunc")),
  1.7,
  tolerance = 1e-10
)

expect_equal(
  as.numeric(family_oichao$variance(eta_oichao, type = "trunc")),
  (1.7 / (3 + 1.7)) * (3 / (3 + 1.7)),
  tolerance = 1e-10
)

expect_equal(
  as.numeric(family_oichao$densityFunction(2, eta_oichao, type = "trunc")),
  dpois(2, 1.7) / (1 - dpois(0, 1.7)),
  tolerance = 1e-10
)

expect_equal(
  as.numeric(family_oichao$densityFunction(3, eta_oichao, type = "nontrunc")),
  dpois(3, 1.7),
  tolerance = 1e-10
)

oichao_dev_y3 <- as.numeric(family_oichao$devResids(
  y = 3, eta = matrix(family_oichao$links[[1]](100)), wt = 1
))
oichao_dev_y2 <- as.numeric(family_oichao$devResids(
  y = 2, eta = matrix(family_oichao$links[[1]](0.01)), wt = 1
))
expect_true(is.finite(oichao_dev_y3))
expect_true(is.finite(oichao_dev_y2))
expect_true(oichao_dev_y3 > 0)
expect_true(oichao_dev_y2 < 0)
expect_true(abs(oichao_dev_y3) < 1)
expect_true(abs(oichao_dev_y2) < 1)

expect_equal(
  as.numeric(family_oichao$devResids(y = 1, eta = eta_oichao, wt = 1)),
  0,
  tolerance = 1e-10
)
expect_equal(
  as.numeric(family_oichao$devResids(y = 4, eta = eta_oichao, wt = 1)),
  0,
  tolerance = 1e-10
)

expect_silent(
  population_oichao_theory <- popSizeEst(fit_oichao_theory)
)
expect_true(is.finite(population_oichao_theory$variance))
expect_true(population_oichao_theory$variance > 0)
expect_true(
  population_oichao_theory$confidenceInterval[1, "lowerBound"] <=
    population_oichao_theory$pointEstimate
)
expect_true(
  population_oichao_theory$confidenceInterval[1, "upperBound"] >=
    population_oichao_theory$pointEstimate
)

expect_error(oichao(lambdaLink = "log"), pattern = "logthird")
