% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.singleRStaticCountData}
\alias{summary.singleRStaticCountData}
\title{Summary statistics for model of singleRStaticCountData class.}
\usage{
\method{summary}{singleRStaticCountData}(
  object,
  test = c("t", "z"),
  resType = "pearson",
  correlation = FALSE,
  confint = FALSE,
  cov,
  popSizeEst,
  ...
)
}
\arguments{
\item{object}{object of singleRStaticCountData class.}

\item{test}{type of test for significance of parameters \code{"t"} for t-test
and \code{"z"} for normal approximation of students t distribution, by
default \code{"z"} is used if there are more than 30 degrees of freedom
and \code{"t"} is used in other cases.}

\item{resType}{type of residuals to summarize any value that is allowed in
\code{residuals.singleRStaticCountData} except for \code{"all"} is allowed. By default
pearson residuals are used.}

\item{correlation}{logical value indicating whether correlation matrix should
be computed from covariance matrix by default \code{FALSE}.}

\item{confint}{logical value indicating whether confidence intervals for
regression parameters should be constructed. By default \code{FALSE}.}

\item{cov}{covariance matrix corresponding to regression parameters.
It is possible to give \code{cov} argument as a function of \code{object}.
If not specified it will be constructed using \code{vcov.singleRStaticCountData} method.
(i.e using Cramer-Rao lower bound)}

\item{popSizeEst}{a \code{popSizeEstResults} class object.
If not specified population size estimation results will be drawn from
\code{object}. If any post-hoc procedures, such as sandwich covariance matrix
estimation or bias reduction, were taken it is possible to include them in
population size estimation results by calling \code{redoPopEstimation}.}

\item{...}{additional optional arguments passed to the following functions:
\itemize{
\item \code{vcov.singleRStaticCountData} -- if no \code{cov} argument was provided.
\item \code{cov} -- if \code{cov} parameter specified at call was a function.
\item \code{confint.singleRStaticCountData} -- if \code{confint} parameter was set to \code{TRUE} at function call.
In particular it is possible to set confidence level in \code{...}.
}}
}
\value{
An object of \code{summarysingleRStaticCountData} class containing:
\itemize{
\item \code{call} -- A call which created \code{object}.
\item \code{coefficients} -- A dataframe with estimated regression coefficients
and their summary statistics such as standard error Wald test statistic and
p value for Wald test.
\item \code{residuals} -- A vector of residuals of type specified at call.
\item \code{aic} -- Akaike's information criterion.
\item \code{bic} -- Bayesian (Schwarz's) information criterion.
\item \code{iter} -- Number of iterations taken in fitting regression.
\item \code{logL} -- Logarithm of likelihood function evaluated at coefficients.
\item \code{deviance} -- Residual deviance.
\item \code{populationSize} -- Object with population size estimation results.
\item \code{dfResidual} -- Residual degrees of freedom.
\item \code{sizeObserved} -- Size of observed population.
\item \code{correlation} -- Correlation matrix if \code{correlation} parameter was set to \code{TRUE}
\item \code{test} -- Type of statistical test performed.
\item \code{model} -- Family class object specified in call for \code{object}.
\item \code{skew} -- If bootstrap sample was saved contains estimate of skewness.
}
}
\description{
A \code{summary} method for \code{singleRStaticCountData} class
}
\details{
\loadmathjax

Works
analogically to \code{summary.glm} but includes population size estimation
results. If any additional statistics, such as confidence intervals for
coefficients or coefficient correlation, are specified they will be printed.
}
\seealso{
\code{\link[=redoPopEstimation]{redoPopEstimation()}} \code{\link[stats:summary.glm]{stats::summary.glm()}}
}
