developer_tests_enabled <- function() {
  isTRUE(tolower(Sys.getenv(
    "TEST_SINGLERCAPTURE_MULTICORE_DEVELOPER",
    unset = "FALSE"
  )) == "true")
}

load_inflated_test_data <- function() {
  list(
    N = 500,
    x1 = rep(0:3, c(99, 228, 137, 36)),
    test_inflated = read.csv("test_inflated.csv")
  )
}

load_zerotruncated_test_data <- function() {
  set.seed(123)

  test_zerotruncated <- read.csv("test_zerotruncated.csv")

  list(
    eps = if (capabilities("long.double")) sqrt(.Machine$double.eps) else 0.01,
    N = 1000,
    test_zerotruncated = test_zerotruncated,
    x1 = test_zerotruncated$x1,
    x2 = test_zerotruncated$x2,
    x3 = test_zerotruncated$x3
  )
}
