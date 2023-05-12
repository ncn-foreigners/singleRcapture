# internal functions

# logit link
singleRinternallogitLink <- function(x,
                                     inverse = FALSE,
                                     deriv = 0) {
  deriv <- deriv + 1
  if (isFALSE(inverse)) {
    res <- switch(
      deriv,
      log(x / (1 - x)),# link
      1 / (x * (1 - x)), # first derivative of link
      (2 * x - 1) / ((x * (1 - x)) ^ 2), # second derivative of link
      -2 * (3 * (x ^ 2) - 3 * x + 1) / ((x * (1 - x)) ^ 3) # third derivative of link
    )
  } else {
    res <- switch(
      deriv,
      1 / (1 + exp(-x)), # inverse link
      exp(x) / ((1 + exp(x)) ^ 2), # first derivative of inverse link
      -exp(x) * (exp(x) - 1) / ((1 + exp(x)) ^ 3), # second derivative of inverse link
      exp(x) * (exp(2 * x) - 4 * exp(x) + 1) / ((1 + exp(x)) ^ 4) # third derivative of inverse link
    )
  }
  
  res
}

# log link
singleRinternallogLink <- function(x,
                                  inverse = FALSE,
                                  deriv = 0) {
  deriv <- deriv + 1
  if (isFALSE(inverse)) {
    res <- switch(
      deriv,
      log(x),# link
      1 / x, # first derivative of link
      -1 / (x ^ 2), # second derivative of link
      2 / (x ^ 3) # third derivative of link
    )
  } else {
    res <- switch(
      deriv,
      exp(x), # inverse link
      exp(x), # first derivative of inverse link
      exp(x), # second derivative of inverse link
      exp(x) # third derivative of inverse link
    )
  }
  
  res
}

# half log link for chao
singleRinternalloghalfLink <- function(x,
                                   inverse = FALSE,
                                   deriv = 0) {
  deriv <- deriv + 1
  if (isFALSE(inverse)) {
    res <- switch(
      deriv,
      log(x / 2),# link
      1 / x, # first derivative of link
      -1 / (x ^ 2), # second derivative of link
      2 / (x ^ 3) # third derivative of link
    )
  } else {
    res <- switch(
      deriv,
      2 * exp(x), # inverse link
      2 * exp(x), # first derivative of inverse link
      2 * exp(x), # second derivative of inverse link
      2 * exp(x) # third derivative of inverse link
    )
  }
  
  res
}

# cloglog

singleRinternalcloglogLink <- function(x,
                                     inverse = FALSE,
                                     deriv = 0) {
  deriv <- deriv + 1
  if (isFALSE(inverse)) {
    res <- switch(
      deriv,
      log(-log(1 - x)),# link
      -1 / ((1 - x) * log(1 - x)), # first derivative of link
      -(1 + log(1 - x)) / ((x - 1) ^ 2 * log(1 - x) ^ 2), # second derivative of link
      (2*log(1 - x) ^ 2 + 3 * log(1 - x) + 2) / (log(1 - x) ^ 3 * (x - 1) ^ 3) # third derivative of link
    )
  } else {
    res <- switch(
      deriv,
      1 - exp(-exp(x)), # inverse link
      exp(x - exp(x)), # first derivative of inverse link
      (1 - exp(x)) * exp(x - exp(x)), # second derivative of inverse link
      (exp(2 * x) - 3 * exp(x) + 1) * exp(x - exp(x)) # third derivative of inverse link
    )
  }
  
  res
}

# probit
#' @importFrom stats dnorm
#' @importFrom stats pnorm
singleRinternalprobitLink <- function(x,
                                      inverse = FALSE,
                                      deriv = 0) {
  deriv <- deriv + 1
  if (isFALSE(inverse)) {
    res <- switch(deriv,
      qnorm(x),# link
      1/dnorm(qnorm(x)), # first derivative of link
      qnorm(x) / (dnorm(qnorm(x))) ^ 2, # second derivative of link
      (1 + 2 * qnorm(x) ^ 2) / dnorm(qnorm(x)) ^ 3 # third derivative of link
    )
  } else {
    res <- switch(deriv,
      pnorm(x), # inverse link
      dnorm(x), # first derivative of inverse link
      -(x * exp(-x ^ 2 / 2)) / (2 * pi) ^ .5, # second derivative of inverse link
      ((x ^ 2 - 1) * exp(-x ^ 2 / 2)) / (2 * pi) ^ .5 # third derivative of inverse link
    )
  }
  
  res
}


# neglog link
singleRinternalneglogLink <- function(x,
                                     inverse = FALSE,
                                     deriv = 0) {
  deriv <- deriv + 1
  if (isFALSE(inverse)) {
    res <- switch(
      deriv,
      -log(x),# link
      -1 / x, # first derivative of link
      1 / x ^ 2, # second derivative of link
      -2 / x ^ 3 # third derivative of link
    )
  } else {
    res <- switch(
      deriv,
      exp(-x), # inverse link
      -exp(-x), # first derivative of inverse link
      exp(-x), # second derivative of inverse link
      -exp(-x) # third derivative of inverse link
    )
  }
  
  res
}