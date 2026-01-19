# Generating data in singleRcapture

An S3 method for
[`stats::simulate`](https://rdrr.io/r/stats/simulate.html) to handle
`singleRStaticCountData` and `singleRfamily` classes.

## Usage

``` r
# S3 method for class 'singleRStaticCountData'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'singleRfamily'
simulate(object, nsim, seed = NULL, eta, truncated = FALSE, ...)
```

## Arguments

- object:

  an object representing a fitted model.

- nsim:

  a numeric scalar specifying:

  - number of response vectors to simulate in
    `simulate.singleRStaticCountData`, defaults to `1L`.

  - number of units to draw in `simulate.singleRfamily`, defaults to
    `NROW(eta)`.

- seed:

  an object specifying if and how the random number generator should be
  initialized (‘seeded’).

- ...:

  additional optional arguments.

- eta:

  a matrix of linear predictors

- truncated:

  logical value indicating whether to sample from truncated or full
  distribution.

## Value

a `data.frame` with `n` rows and `nsim` columns.

## See also

[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html)
[`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md)

## Author

Maciej Beręsewicz, Piotr Chlebicki

## Examples

``` r
N <- 10000
###gender <- rbinom(N, 1, 0.2)
gender <- rep(0:1, c(8042, 1958))
eta <- -1 + 0.5*gender
counts <- simulate(ztpoisson(), eta = cbind(eta), seed = 1)
df <- data.frame(gender, eta, counts)
df2 <- subset(df, counts > 0)
### check coverage with summary
mod1 <-  estimatePopsize(
  formula       = counts ~ 1 + gender, 
  data          = df2, 
  model         = ztpoisson, 
  controlMethod = list(silent = TRUE)
)
mod1_sims <- simulate(mod1, nsim=10, seed = 1)
colMeans(mod1_sims)
#>    sim_1    sim_2    sim_3    sim_4    sim_5    sim_6    sim_7    sim_8 
#> 1.241014 1.259281 1.239246 1.244255 1.226871 1.237773 1.239540 1.234532 
#>    sim_9   sim_10 
#> 1.230701 1.244844 
mean(df2$counts)
#> [1] 1.240424
```
