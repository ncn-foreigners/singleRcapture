# Diagnostic plots for regression and population size estimation

Simple diagnostic plots for `singleRStaticCountData` class objects.

## Usage

``` r
# S3 method for class 'singleRStaticCountData'
plot(
  x,
  plotType = c("qq", "marginal", "fitresid", "bootHist", "rootogram", "dfpopContr",
    "dfpopBox", "scaleLoc", "cooks", "hatplot", "strata"),
  confIntStrata = c("normal", "logNormal"),
  histKernels = TRUE,
  dfpop,
  ...
)
```

## Arguments

- x:

  object of `singleRStaticCountData` class.

- plotType:

  character parameter (default `"qq"`) specifying type of plot to be
  made. The following list presents and briefly explains possible type
  of plots:

  - `qq` – The quantile-quantile plot for pearson residuals (or
    standardized pearson residuals if these are available for the model)
    i.e. empirical quantiles from residuals are plotted against
    theoretical quantiles from standard distribution.

  - `marginal` – A plot made by `matplot` with fitted and observed
    marginal frequencies with labels.

  - `fitresid` – Plot of fitted linear predictors against (standardized)
    pearson residuals.

  - `bootHist` – Simple histogram for statistics obtained from
    bootstrapping (if one was performed and the statistics were saved).

  - `rootogram` – Rootogram, for full explanation see: Kleiber and
    Zeileis Visualizing Count Data Regressions Using Rootograms (2016),
    in short it is a `barplot` where height is the square root of
    observed marginal frequencies adjusted by difference between square
    root of observed and fitted marginal frequencies connected by line
    representing fitted marginal frequencies. The less of a difference
    there is between the 0 line and beginning of a bar the more accurate
    fitt was produced by the model.

  - `dfpopContr` – Plot of `dfpopsize` against unit contribution. On the
    plot is y = x line i.e. what deletion effect would be if removing
    the unit from the model didn't effect regression coefficients. The
    further away the observation is from this line the more influential
    it is.

  - `dfpopBox` – Boxplot of `dfpopsize` for getting the general idea
    about the distribution of the "influence" of each unit on population
    size estimate.

  - `scaleLoc` – The scale - location plot i.e. square root of absolute
    values of (standardized) pearson residuals against linear predictors
    for each column of linear predictors.

  - `cooks` – Plot of cooks distance for detecting influential
    observations.

  - `hatplot` – Plot of hat values for each linear predictor for
    detecting influential observations.

  - `strata` – Plot of confidence intervals and point estimates for
    strata provided in `...` argument

- confIntStrata:

  confidence interval type to use for strata plot. Currently supported
  values are `"normal"` and `"logNormal"`.

- histKernels:

  logical value indicating whether to add density lines to histogram.

- dfpop:

  TODO

- ...:

  additional optional arguments passed to the following functions:

  - For `plotType = "bootHist"`

    - [`graphics::hist`](https://rdrr.io/r/graphics/hist.html) – with
      `x, main, xlab, ylab` parameters fixed.

  - For `plotType = "rootogram"`

    - [`graphics::barplot`](https://rdrr.io/r/graphics/barplot.html) –
      with `height, offset, ylab, xlab, ylim` parameters fixed.

    - [`graphics::lines`](https://rdrr.io/r/graphics/lines.html) – with
      `x, y, pch, type, lwd, col` parameters fixed.

  - For `plotType = "dfpopContr"`

    - `dfpopsize` – with `model, observedPop` parameters fixed.

    - `plot.default` – with `x, y, xlab, main` parameters fixed.

  - For `plotType = "dfpopBox"`

    - `dfpopsize` – with `model, observedPop` parameters fixed.

    - [`graphics::boxplot`](https://rdrr.io/r/graphics/boxplot.html) –
      with `x, ylab, main` parameters fixed.

  - For `plotType = "scaleLoc"`

    - `plot.default` – with `x, y, xlab, ylab, main, sub` parameters
      fixed.

  - For `plotType = "fitresid"`

    - `plot.default` – with `x, y, xlab, ylab, main, sub` parameters
      fixed.

  - For `plotType = "cooks"`

    - `plot.default` – with `x, xlab, ylab, main` parameters fixed.

  - For `plotType = "hatplot"`

    - `hatvalues.singleRStaticCountData`

    - `plot.default` – with `x, xlab, ylab, main` parameters fixed.

  - For `plotType = "strata"`

    - `stratifyPopsize.singleRStaticCountData`

## Value

No return value only the plot being made.

## See also

[`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md)
[`dfpopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/regDiagSingleR.md)
[`marginalFreq()`](https://ncn-foreigners.github.io/singleRcapture/reference/marginalFreq.md)
[`stats::plot.lm()`](https://rdrr.io/r/stats/plot.lm.html)
[`stats::cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)
[`hatvalues.singleRStaticCountData()`](https://ncn-foreigners.github.io/singleRcapture/reference/regDiagSingleR.md)

## Author

Piotr Chlebicki
