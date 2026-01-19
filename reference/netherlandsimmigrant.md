# Immigration data in the Netherlands

This data set contains information about immigrants in four cities
(Amsterdam, Rotterdam, The Hague and Utrecht) in Netherlands that have
been staying in the country illegally in 1995 and have appeared in
police records that year.

## Usage

``` r
data("netherlandsimmigrant")
```

## Format

Data frame with 1,880 rows and 5 columns.

- `capture`:

  Number of times a person has been captured by police.

- `gender`:

  Factor describing gender of the apprehended person.

- `age`:

  Factor describing age of apprehended person. Either bellow or above 40
  years old.

- `reason`:

  Factor describing reason for being apprehended by police either
  illegal stay in Netherlands or other reasons.

- `nation`:

  Factor with nation of origin of the captured person. There are 6
  levels of this variable:
  `"American and Australia", "Asia", "North Africa", "Rest of Africa", "Surinam", "Turkey"`.

## References

This data set and its description was provided in publication: van Der
Heijden, P. G., Bustami, R., Cruyff, M. J., Engbersen, G., and Van
Houwelingen, H. C. (2003). Point and interval estimation of the
population size using the truncated Poisson regression model.
*Statistical Modelling*, 3(4), 305-322. doi:10.1191/1471082X03st057oa
