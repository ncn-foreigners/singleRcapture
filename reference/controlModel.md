# Control parameters specific to some models

`controlModel` constructs a list with all necessary control parameters
in `estimatePopsize` that are either specific to selected model or do
not fit anywhere else.

Specifying additional formulas should be done by using only right hand
side of the formula also for now all variables from additional formulas
should also be included in the "main" formula.

## Usage

``` r
controlModel(
  weightsAsCounts = FALSE,
  omegaFormula = ~1,
  alphaFormula = ~1,
  piFormula = ~1
)
```

## Arguments

- weightsAsCounts:

  a boolean value indicating whether to treat `weights` argument as
  number of occurrences for each row in the `data` and adjust necessary
  methods and functionalities, like adjustments in bootstrap or
  decreasing weights in `dfbeta` instead or deleting rows from data, to
  accommodate this form of model specification.

- omegaFormula:

  a formula for inflation parameter in one inflated zero truncated and
  zero truncated one inflated models.

- alphaFormula:

  a formula for dispersion parameter in negative binomial based models.

- piFormula:

  a formula for probability parameter in pseudo hurdle zero truncated
  and zero truncated pseudo hurdle models.

## Value

A list with selected parameters, it is also possible to call list
directly.

## See also

[`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md)
[`controlMethod()`](https://ncn-foreigners.github.io/singleRcapture/reference/controlMethod.md)
[`controlPopVar()`](https://ncn-foreigners.github.io/singleRcapture/reference/controlPopVar.md)
[`singleRmodels()`](https://ncn-foreigners.github.io/singleRcapture/reference/singleRmodels.md)

## Author

Piotr Chlebicki, Maciej Beręsewicz
