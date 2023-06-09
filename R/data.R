#' @title British farm submissions data
#' @description Data on British animal farms submissions to AHVLA. British farms 
#' are able to submit samples to AHVLA if cause of death for an animal
#' cannot be determined and private veterinary surgeon decides to submit them, 
#' unless there is notifiable disease suspected then such a submission is not 
#' required. 
#' 
#' This data set contains information about such farms.
#' All submissions from farms are included in this data frame
#' not only carcasses but also blood samples etc.
#' @docType data
#' @format Data frame with 12,036 rows and 4 columns.
#' \describe{
#'   \item{\code{TOTAL_SUB}}{Number of submissions of animal material.}
#'   \item{\code{log_size}}{Numerical value equal to logarithm of size of farm.}
#'   \item{\code{log_distance}}{Numerical value equal to logarithm of distance to nearest AHVLA center.}
#'   \item{\code{C_TYPE}}{Factor describing type of activity on farm that animals are 
#'   used for. Either \code{Dairy} or \code{Beef}}
#' }
#' @name farmsubmission
#' @references This data set and its description was provided in publication:
#' B&ouml;hning, D., Vidal Diez, A., Lerdsuwansri, R., Viwatwongkasem, C., and Arnold, M. (2013). "A generalization of Chao's estimator for covariate information". *Biometrics*, 69(4), 1033-1042. doi:10.1111/biom.12082
#' @usage data("farmsubmission")
"farmsubmission"

#' @title Data on immigration in Netherlands
#' @description  This data set contains information about immigrants in four
#' cities (Amsterdam, Rotterdam, The Hague and Utrecht) in Netherlands that have
#' been staying in the country illegally in 1995 and have appeared in police
#' records that year.
#' @docType data
#' @format Data frame with 1,880 rows and 5 columns.
#' \describe{
#'   \item{\code{capture}}{Number of times a person has been captured by police.}
#'   \item{\code{gender}}{Factor describing gender of the apprehended person.}
#'   \item{\code{age}}{Factor describing age of apprehended person. 
#'   Either bellow or above 40 years old.}
#'   \item{\code{reason}}{Factor describing reason for being apprehended by 
#'   police either illegal stay in Netherlands or other reasons.}
#'   \item{\code{nation}}{Factor with nation of origin of the captured person. 
#'   There are 6 levels of this variable: \code{"American and Australia",
#'   "Asia", "North Africa", "Rest of Africa", "Surinam", "Turkey"}.}
#' }
#' @name netherlandsimmigrant
#' @references This data set and its description was provided in publication:
#' van Der Heijden, P. G., Bustami, R., Cruyff, M. J., Engbersen, G., and Van Houwelingen, H. C. (2003). Point and interval estimation of the population size using the truncated Poisson regression model. *Statistical Modelling*, 3(4), 305-322. doi:10.1191/1471082X03st057oa
#' @usage data("netherlandsimmigrant")
"netherlandsimmigrant"

#' @title British farm carcass submissions data
#' @description Data on British animal farms submissions to AHVLA. British farms 
#' are able to submit samples to AHVLA if cause of death for an animal
#' cannot be determined and private veterinary surgeon decides to submit them, 
#' unless there is notifiable disease suspected then such a submission is not 
#' required. 
#' 
#' This data set contains information about such farms.
#' Only submissions that are included in this data frame
#' are submissions of carcasses i.e. submissions of blood samples etc. 
#' are excluded.
#' @docType data
#' @format Data frame with 1,858 rows and 4 columns.
#' \describe{
#'   \item{\code{TOTAL_SUB}}{Number of submissions of animal carcasses.}
#'   \item{\code{log_size}}{Numerical value equal to logarithm of size of farm.}
#'   \item{\code{log_distance}}{Numerical value equal to logarithm of distance to nearest AHVLA center.}
#'   \item{\code{C_TYPE}}{Factor describing type of activity on farm that animals are 
#'   used for. Either \code{Dairy} or \code{Beef}}
#' }
#' @name carcassubmission
#' @references This data set and its description was provided in publication:
#' B&ouml;hning, D., Vidal Diez, A., Lerdsuwansri, R., Viwatwongkasem, C., and Arnold, M. (2013). "A generalization of Chao's estimator for covariate information". *Biometrics*, 69(4), 1033-1042. doi:10.1111/biom.12082
#' @usage data("carcassubmission")
"carcassubmission"
