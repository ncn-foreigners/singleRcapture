# Observed and fitted marginal frequencies

A function that given a fitted `singleR` class object computed marginal
frequencies by as sum of probability density functions for each unit in
data at each point i.e. kth element of marginal frequency table is given
by \\\sum\_{j=1}^{N\_{obs}}\mathbb{P}(Y\_{j}=k\|\eta\_{j})\\. For k=0
only (if specified at call) they are computed as \\\hat{N}-N\_{obs}\\
because \\\boldsymbol{f}\_{0}\\ is assumed to the unobserved part of the
studied population.

These frequencies are useful in diagnostics for count data regression,
such as assessment of fit.

## Usage

``` r
marginalFreq(
  object,
  includeones = TRUE,
  includezeros = TRUE,
  onecount = NULL,
  range,
  ...
)
```

## Arguments

- object:

  object of `singleR` class.

- includeones:

  logical value indicating whether to include the estimated number of
  zero counts.

- includezeros:

  logical value indicating whether to include one counts in the zero-one
  truncated models.

- onecount:

  a numeric value indicating number of one counts if null `trcount` from
  object will be assumed to be a number one counts.

- range:

  optional argument specifying range of selected Y values.

- ...:

  currently does nothing.

## Value

A list with observed name of the fitted model family degrees of freedom
and observed and fitted marginal frequencies.

## See also

[`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md)
– where example of usage is provided

## Author

Piotr Chlebicki
