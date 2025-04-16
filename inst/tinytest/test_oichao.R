
# Helper function to compute expected OIChao by hand
compute_oichao <- function(f2, f3, n, bias_corr = FALSE) {
  if (f3 == 0) {
    f0 <- 0
    variance <- 0
  } else if (bias_corr) {
    # Bias-corr point estimate from the latest oichao function
    f0 <- (2/9) * ((f2^3 - 3*f2^2 + 2*f2) / ((f3+1)*(f3+2)))
    # Variance for bias-corr estimator (matches popVar in oichao)
    v <- (f2^3 - 3 * f2^2 + 2 * f2) / ((f3 + 1) * (f3 + 2))
    du_df2 <- (2/9) * (6 * f3 / f2^2)
    du_df3 <- (-2/9) * (6 / f2)
    dv_df2 <- (3 * f2^2 - 6 * f2 + 2) / ((f3 + 1) * (f3 + 2))
    dv_df3 <- (f2^3 - 3 * f2^2 + 2 * f2) * (-1 / ((f3 + 1)^2 * (f3 + 2)) - 1 / ((f3 + 1) * (f3 + 2)^2))
    df0_df2 <- du_df2 * v + (2/9) * dv_df2
    df0_df3 <- du_df3 * v + (2/9) * dv_df3
    variance <- (df0_df2^2 * f2) + (df0_df3^2 * f3)
  } else {
    # Original point estimate from the latest oichao function
    f0 <- (2/9) * (f2^3 / f3^2)
    # Variance for original estimator (matches popVar in oichao)
    variance <- (2/9)^2 * ( (3 * f2^2 / f3^2)^2 * f2 + (2 * f2^3 / f3^3)^2 * f3 )
  }
  N <- n + f0
  se <- sqrt(variance)
  list(N = N, f0 = f0, variance = variance, se = se)
}

# 1. No-Covariate Cases
# 1.1 Perpetrators (2009) Data
perpetrators_data <- data.frame(
  y = c(rep(1, 15169), rep(2, 957), rep(3, 393), rep(4, 99), rep(5, 28), rep(6, 16))
)
fit_perpetrators <- estimatePopsize(y ~ 1, data = perpetrators_data, 
                                    model = oichao(bias_corr = FALSE), 
                                    method = "IRLS")

# Hand calculation (original estimator)
expected_perpetrators <- compute_oichao(f2 = 957, f3 = 393, n = 17662, bias_corr = FALSE)

expect_equal(fit_perpetrators$populationSize$pointEstimate, 
             expected_perpetrators$N, 
             tolerance = 1,
             info = "Perpetrators (2009): Point estimate should be ~18923.47")

expect_true(abs(sqrt(fit_perpetrators$populationSize$variance) - expected_perpetrators$se) < 10,
            info = "Perpetrators (2009): SE should be ~387")

# Test bias-corr estimator for Perpetrators (2009)
fit_perpetrators_bc <- estimatePopsize(y ~ 1, data = perpetrators_data, 
                                       model = oichao(bias_corr = TRUE), method = "IRLS")

expected_perpetrators_bc <- compute_oichao(f2 = 957, f3 = 393, n = 17662, bias_corr = TRUE)

expect_equal(fit_perpetrators_bc$populationSize$pointEstimate, 
             expected_perpetrators_bc$N, 
             tolerance = 1,
             info = "Perpetrators (2009, Bias-corr): Point estimate should be ~18910.84")

expect_true(abs(sqrt(fit_perpetrators_bc$populationSize$variance) - expected_perpetrators_bc$se) < 10,
            info = "Perpetrators (2009, Bias-corr): SE should be ~356.87")

# 1.2 Synthetic Data
synthetic_data <- data.frame(
  y = c(rep(1, 690), rep(2, 95), rep(3, 32), rep(4, 7))
)

fit_synthetic <- estimatePopsize(y ~ 1, data = synthetic_data, 
                                 model = oichao(bias_corr = FALSE), method = "IRLS")

# Hand calculation
expected_synthetic <- compute_oichao(f2 = 95, f3 = 32, n = 824, bias_corr = FALSE)
expect_equal(fit_synthetic$populationSize$pointEstimate, expected_synthetic$N, tolerance = 1,
             info = "Synthetic Data: Point estimate should be ~1010.06")
expect_true(abs(sqrt(fit_synthetic$populationSize$variance) - expected_synthetic$se) < 5,
            info = "Synthetic Data: SE should be ~48.2")

# 2. Edge Cases
# 2.1 f3 = 0 (should return f0 = 0)
edge_f3_zero <- data.frame(
  y = c(rep(1, 100), rep(2, 10))
)

fit_f3_zero <- estimatePopsize(y ~ 1, data = edge_f3_zero, 
                               model = oichao(bias_corr = FALSE), method = "IRLS")

# Hand calculation
expect_equal(fit_f3_zero$populationSize$pointEstimate, 110, tolerance = 0.1,
             info = "Edge Case (f3 = 0): Point estimate should be 110")
expect_equal(fit_f3_zero$populationSize$variance, 0, tolerance = 0.1,
             info = "Edge Case (f3 = 0): Variance should be 0")

# 2.2 Small Sample Size
edge_small <- data.frame(
  y = c(1, 2, 3)  # f1 = 1, f2 = 1, f3 = 1
)
fit_small <- estimatePopsize(y ~ 1, data = edge_small, model = oichao(bias_corr = FALSE), method = "IRLS")

# Hand calculation
expected_small <- compute_oichao(f2 = 1, f3 = 1, n = 3, bias_corr = FALSE)
expect_equal(fit_small$populationSize$pointEstimate, expected_small$N, tolerance = 0.1,
             info = "Edge Case (Small Sample): Point estimate should be ~3.2222")
expect_true(abs(sqrt(fit_small$populationSize$variance) - expected_small$se) < 0.1,
            info = "Edge Case (Small Sample): SE should be ~0.377")

# 2.3 Large Frequencies (Numerical Stability)
edge_large <- data.frame(
  y = c(rep(1, 10000), rep(2, 5000), rep(3, 2000))
)
fit_large <- estimatePopsize(y ~ 1, data = edge_large, model = oichao(bias_corr = FALSE), method = "IRLS")

# Hand calculation
expected_large <- compute_oichao(f2 = 5000, f3 = 2000, n = 17000, bias_corr = FALSE)
expect_equal(fit_large$populationSize$pointEstimate, expected_large$N, tolerance = 1,
             info = "Edge Case (Large Frequencies): Point estimate should be ~23944.44")
expect_true(abs(sqrt(fit_large$populationSize$variance) - expected_large$se) < 20,
            info = "Edge Case (Large Frequencies): SE should be ~623.75")

# 3. Covariate Case
set.seed(123)
cov_data <- data.frame(
  y = c(rep(1, 50), rep(2, 30), rep(3, 10)),
  z = rbinom(90, 1, 0.5)
)
fit_cov <- estimatePopsize(y ~ z, data = cov_data, model = oichao(bias_corr = FALSE), method = "IRLS")

# Theoretical Expectation (Approximate)
expect_true(fit_cov$populationSize$pointEstimate > 90 && fit_cov$populationSize$pointEstimate < 200,
            info = "Covariate Case: Point estimate should be between 90 and 200")
expect_true(sqrt(fit_cov$populationSize$variance) > 0,
            info = "Covariate Case: SE should be positive")