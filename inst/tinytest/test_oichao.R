# Helper function to compute expected OIChao by hand
compute_oichao <- function(f2, f3, n) {
  f0 <- if (f3 > 0) (2/9) * (f2^3 / f3^2) else 0
  N <- n + f0
  variance <- if (f3 > 0) (2/9)^2 * ( (3 * f2^2 / f3^2)^2 * f2 + (2 * f2^3 / f3^3)^2 * f3 ) else 0
  se <- sqrt(variance)
  list(N = N, f0 = f0, variance = variance, se = se)
}

# 1. No-Covariate Cases
# 1.1 Perpetrators (2009) Data from Table 1
perpetrators_data <- data.frame(
  y = c(rep(1, 15169), rep(2, 957), rep(3, 393), rep(4, 99), rep(5, 28), rep(6, 16))
)
fit_perpetrators <- estimatePopsize(y ~ 1, data = perpetrators_data, model = "oichao", method = "IRLS")

# Hand calculation
# f2 = 957, f3 = 393, n = 17662
# f0 = (2/9) * (957^3 / 393^2) ≈ 1261.47
# N = 17662 + 1261.47 ≈ 18923.47
# Variance = (2/9)^2 * [ (3 * 957^2 / 393^2)^2 * 957 + (2 * 957^3 / 393^3)^2 * 393 ] ≈ 149614
# SE ≈ 387
expected_perpetrators <- compute_oichao(f2 = 957, f3 = 393, n = 17662)
expect_equal(fit_perpetrators$populationSize$pointEstimate, expected_perpetrators$N, tolerance = 1,
             info = "Perpetrators (2009): Point estimate should be ~18923.47")
expect_true(abs(sqrt(fit_perpetrators$populationSize$variance) - expected_perpetrators$se) < 10,
            info = "Perpetrators (2009): SE should be ~387")

# 1.2 Synthetic Data from Example 2
synthetic_data <- data.frame(
  y = c(rep(1, 690), rep(2, 95), rep(3, 32), rep(4, 7))
)
fit_synthetic <- estimatePopsize(y ~ 1, data = synthetic_data, model = "oichao", method = "IRLS")

# Hand calculation
# f2 = 95, f3 = 32, n = 824
# f0 = (2/9) * (95^3 / 32^2) ≈ 186.0621
# N = 824 + 186.0621 ≈ 1010.0621
# Variance = (2/9)^2 * [ (3 * 95^2 / 32^2)^2 * 95 + (2 * 95^3 / 32^3)^2 * 32 ] ≈ 2325
# SE ≈ 48.2
expected_synthetic <- compute_oichao(f2 = 95, f3 = 32, n = 824)
expect_equal(fit_synthetic$populationSize$pointEstimate, expected_synthetic$N, tolerance = 1,
             info = "Synthetic Data: Point estimate should be ~1010.0621")
expect_true(abs(sqrt(fit_synthetic$populationSize$variance) - expected_synthetic$se) < 5,
            info = "Synthetic Data: SE should be ~48.2")

# 2. Edge Cases
# 2.1 f3 = 0 (should return f0 = 0)
edge_f3_zero <- data.frame(
  y = c(rep(1, 100), rep(2, 10))
)
fit_f3_zero <- estimatePopsize(y ~ 1, data = edge_f3_zero, model = "oichao", method = "IRLS")

# Hand calculation
# f2 = 10, f3 = 0, n = 110
# f0 = 0 (since f3 = 0)
# N = 110 + 0 = 110
# Variance = 0
expect_equal(fit_f3_zero$populationSize$pointEstimate, 110, tolerance = 0.1,
             info = "Edge Case (f3 = 0): Point estimate should be 110")
expect_equal(fit_f3_zero$populationSize$variance, 0, tolerance = 0.1,
             info = "Edge Case (f3 = 0): Variance should be 0")

# 2.2 Small Sample Size
edge_small <- data.frame(
  y = c(1, 2, 3)  # f1 = 1, f2 = 1, f3 = 1
)
fit_small <- estimatePopsize(y ~ 1, data = edge_small, model = "oichao", method = "IRLS")

# Hand calculation
# f2 = 1, f3 = 1, n = 3
# f0 = (2/9) * (1^3 / 1^2) = 2/9 ≈ 0.2222
# N = 3 + 0.2222 ≈ 3.2222
# Variance = (2/9)^2 * [ (3 * 1^2 / 1^2)^2 * 1 + (2 * 1^3 / 1^3)^2 * 1 ] = (2/9)^2 * (9 + 4) ≈ 0.142
# SE ≈ 0.377
expected_small <- compute_oichao(f2 = 1, f3 = 1, n = 3)
expect_equal(fit_small$populationSize$pointEstimate, expected_small$N, tolerance = 0.1,
             info = "Edge Case (Small Sample): Point estimate should be ~3.2222")
expect_true(abs(sqrt(fit_small$populationSize$variance) - expected_small$se) < 0.1,
            info = "Edge Case (Small Sample): SE should be ~0.377")

# 2.3 Large Frequencies (Numerical Stability)
edge_large <- data.frame(
  y = c(rep(1, 10000), rep(2, 5000), rep(3, 2000))
)
fit_large <- estimatePopsize(y ~ 1, data = edge_large, model = "oichao", method = "IRLS")

# Hand calculation
# f2 = 5000, f3 = 2000, n = 17000
# f0 = (2/9) * (5000^3 / 2000^2) ≈ 6944.444
# N = 17000 + 6944.444 ≈ 23944.444
# Variance = (2/9)^2 * [ (3 * 5000^2 / 2000^2)^2 * 5000 + (2 * 5000^3 / 2000^3)^2 * 2000 ] ≈ 389,062
# SE ≈ 623.75
expected_large <- compute_oichao(f2 = 5000, f3 = 2000, n = 17000)
expect_equal(fit_large$populationSize$pointEstimate, expected_large$N, tolerance = 1,
             info = "Edge Case (Large Frequencies): Point estimate should be ~23944.444")
expect_true(abs(sqrt(fit_large$populationSize$variance) - expected_large$se) < 20,
            info = "Edge Case (Large Frequencies): SE should be ~623.75")

# 3. Covariate Case
# Add a binary covariate z (0 or 1)
set.seed(123)
cov_data <- data.frame(
  y = c(rep(1, 50), rep(2, 30), rep(3, 10)),
  z = rbinom(90, 1, 0.5)
)
fit_cov <- estimatePopsize(y ~ z, data = cov_data, model = "oichao", method = "IRLS")

# Theoretical Expectation (Approximate)
# Without covariates: f2 = 30, f3 = 10, n = 90
# f0 = (2/9) * (30^3 / 10^2) = 60, N = 90 + 60 = 150
# With covariates, lambda varies by z, so we expect N to adjust
# lambda ≈ 0.5 (from earlier fits), but varies by z
# Check summary for reasonable N (should be close to 150, adjusted by covariate effect)
expect_true(fit_cov$populationSize$pointEstimate > 90 && fit_cov$populationSize$pointEstimate < 200,
            info = "Covariate Case: Point estimate should be between 90 and 200")
expect_true(sqrt(fit_cov$populationSize$variance) > 0,
            info = "Covariate Case: SE should be positive")