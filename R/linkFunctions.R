# internal functions

# logit link
singleRinternallogitLink <- function(x,
                                     inverse = FALSE,
                                     deriv = 0) {
  deriv <- deriv + 1
  if (isFALSE(inverse)) {
    res <- switch(
      deriv,
      log(x / (1 - x)),
      1 / (x * (1 - x)),
      (2 * x - 1) / ((x * (1 - x)) ^ 2), 
      -2 * (3 * (x ^ 2) - 3 * x + 1) / ((x * (1 - x)) ^ 3)
    )
  } else {
    res <- switch(
      deriv,
      1 / (1 + exp(-x)),
      exp(x) / ((1 + exp(x)) ^ 2),
      -exp(x) * (exp(x) - 1) / ((1 + exp(x)) ^ 3),
      exp(x) * (exp(2 * x) - 4 * exp(x) + 1) / ((1 + exp(x)) ^ 4)
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
      log(x),
      1 / x, 
      -1 / (x ^ 2), 
      2 / (x ^ 3)
    )
  } else {
    res <- exp(x)
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
      log(x / 2),
      1 / x, 
      -1 / (x ^ 2),
      2 / (x ^ 3) 
    )
  } else {
    res <- 2 * exp(x)
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
      log(-log(1 - x)),
      -1 / ((1 - x) * log(1 - x)),
      -(1 + log(1 - x)) / ((x - 1) ^ 2 * log(1 - x) ^ 2),
      (2*log(1 - x) ^ 2 + 3 * log(1 - x) + 2) / (log(1 - x) ^ 3 * (x - 1) ^ 3)
    )
  } else {
    res <- switch(
      deriv,
      1 - exp(-exp(x)),
      exp(x - exp(x)),
      (1 - exp(x)) * exp(x - exp(x)),
      (exp(2 * x) - 3 * exp(x) + 1) * exp(x - exp(x))
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
      qnorm(x),
      1/dnorm(qnorm(x)),
      qnorm(x) / (dnorm(qnorm(x))) ^ 2,
      (1 + 2 * qnorm(x) ^ 2) / dnorm(qnorm(x)) ^ 3
    )
  } else {
    res <- switch(deriv,
      pnorm(x),
      dnorm(x),
      -(x * exp(-x ^ 2 / 2)) / (2 * pi) ^ .5,
      ((x ^ 2 - 1) * exp(-x ^ 2 / 2)) / (2 * pi) ^ .5
    )
  }
  
  res
}


# neglog link
singleRinternalneglogLink <- function(x,
                                     inverse = FALSE,
                                     deriv = 0) {
  if (isFALSE(inverse)) {
    deriv <- deriv + 1
    res <- switch(
      deriv,
      -log(x),
      -1 / x, 
      1 / x ^ 2, 
      -2 / x ^ 3 
    )
  } else {
    res <- exp(-x) * ((-1) ^ deriv)
  }
  
  res
}