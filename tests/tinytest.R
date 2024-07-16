
if ( requireNamespace("tinytest", quietly=TRUE) ){
  xx <- isTRUE(tolower(Sys.getenv("TEST_SINGLERCAPTURE_MULTICORE_DEVELOPER")) == "true")
  tinytest::test_package("singleRcapture", ncpu = if (xx) 2 else NULL)
}

