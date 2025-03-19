# test_ratioReg.R
library(testthat)

test_that("Testing ratioReg function with small dataset", {
  test_data <- data.frame(observed = c(1, 1, 2, 2, 2, 3, 3, 3, 4, 5))
  
  result <- ratioReg(data = test_data, formula = observed ~ 1)
  
  expect_true(is.list(result))
  expect_true("estimates" %in% names(result))
  expect_true("confidenceInterval" %in% names(result))
  expect_equal(length(result$estimates), 2)
  expect_true(all(!is.na(result$estimates)))
})

test_that("Testing ratioReg function with negative values", {
  test_data <- data.frame(observed = c(-1, 1, 2, 2, 3, 3, 3, 4, 5))
  
  expect_error(ratioReg(data = test_data, formula = observed ~ 1), 
               "Negative values in 'observed' data are not allowed.")
})

test_that("Testing ratioReg function with larger dataset", {
  test_data <- data.frame(observed = sample(1:100, 100, replace = TRUE))
  
  result <- ratioReg(data = test_data, formula = observed ~ 1)
  
  expect_true("estimates" %in% names(result))
  expect_true("confidenceInterval" %in% names(result))
  expect_true(length(result$estimates) == 2)
})

test_that("Testing ratioReg function with zero counts", {
  test_data <- data.frame(observed = c(0, 1, 2, 3, 4))
  
  expect_error(ratioReg(data = test_data, formula = observed ~ 1), 
               "Zero counts in the data are not allowed.")
})
