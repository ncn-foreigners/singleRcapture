# The British farm submissions data

Data on British animal farms submissions to AHVLA. British farms are
able to submit samples to AHVLA if cause of death for an animal cannot
be determined and private veterinary surgeon decides to submit them,
unless there is notifiable disease suspected then such a submission is
not required.

This data set contains information about such farms. All submissions
from farms are included in this data frame not only carcasses but also
blood samples etc.

## Usage

``` r
data("farmsubmission")
```

## Format

Data frame with 12,036 rows and 4 columns.

- `TOTAL_SUB`:

  Number of submissions of animal material.

- `log_size`:

  Numerical value equal to logarithm of size of farm.

- `log_distance`:

  Numerical value equal to logarithm of distance to nearest AHVLA
  center.

- `C_TYPE`:

  Factor describing type of activity on farm that animals are used for.
  Either `Dairy` or `Beef`

## References

This data set and its description was provided in publication: Böhning,
D., Vidal Diez, A., Lerdsuwansri, R., Viwatwongkasem, C., and Arnold, M.
(2013). "A generalization of Chao's estimator for covariate
information". *Biometrics*, 69(4), 1033-1042. doi:10.1111/biom.12082
