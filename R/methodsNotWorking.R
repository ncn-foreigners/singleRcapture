# Add in the future actual methods for this:
#' @importFrom stats add1
#' @method add1 singleRStaticCountData
#' @exportS3Method 
add1.singleRStaticCountData <- function(object, scope, ...) {
  stop("The add1 method for singleRStaticCountData class doesn't work yet.")
}
#' @importFrom stats profile
#' @method profile singleRStaticCountData
#' @exportS3Method 
profile.singleRStaticCountData <- function(fitted, ...) {
  stop("The profile method for singleRStaticCountData class doesn't work yet.")
}
#' @importFrom stats drop1
#' @method drop1 singleRStaticCountData
#' @exportS3Method 
drop1.singleRStaticCountData <- function(object, scope, ...) {
  stop("The drop1 method for singleRStaticCountData class doesn't work yet.")
}
#' @importFrom stats anova
#' @method anova singleRStaticCountData
#' @exportS3Method 
anova.singleRStaticCountData <- function(object, ...) {
  stop("The anova method for singleRStaticCountData class doesn't work yet.")
}
