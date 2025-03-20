test_data <- data.frame(observed = c(1, 1, 2, 2, 2, 3, 3, 3, 4, 5))

result <- ratioReg(data = test_data, formula = observed ~ 1)

expect_true(is.list(result), "Result should be a list")
expect_true("estimates" %in% names(result), "Result should contain 'estimates'")
expect_true("confidenceInterval" %in% names(result), "Result should contain 'confidenceInterval'")
expect_equal(length(result$estimates), 2, tolerance = 1e-6)
expect_true(all(!is.na(result$estimates)), "All estimates should be non-NA")

# Test negative values
test_data_neg <- data.frame(observed = c(-1, 1, 2, 2, 3, 3, 3, 4, 5))

expect_error(
  ratioReg(data = test_data_neg, formula = observed ~ 1), 
  "Negative values in 'observed' data are not allowed.",
  info = "Function should return an error for negative values"
)

# Test larger dataset
test_data_large <- data.frame(observed = sample(1:100, 100, replace = TRUE))

result_large <- ratioReg(data = test_data_large, formula = observed ~ 1)

expect_true("estimates" %in% names(result_large), "Large dataset should contain 'estimates'")
expect_true("confidenceInterval" %in% names(result_large), "Large dataset should contain 'confidenceInterval'")
expect_equal(length(result_large$estimates), 2, tolerance = 1e-6)

# Test zero counts
test_data_zero <- data.frame(observed = c(0, 1, 2, 3, 4))

expect_error(
  ratioReg(data = test_data_zero, formula = observed ~ 1), 
  "Zero counts in the data are not allowed.",
  info = "Function should return an error for zero values"
)
