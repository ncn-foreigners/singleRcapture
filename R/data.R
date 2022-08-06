#' @title All farm submissions
#' @details Data on British animal farms submissions to AHVLA. British farms 
#' are able to submit samples to AHVLA if cause of death of an animal
#' cannot be determined and private veterinary surgeon decides to submit them, 
#' unless there is notifiable disease suspected then it is not required. 
#' This is data on such farms and their submissions. 
#' This data set contains farms that have submitted anything at all 
#' not only carcasses but also blood samples etc.
#' @docType data
#' @format Data frame with 12036 rows and 4 columns.
#' \describe{
#'   \item{TOTAL_SUB}{Number of submissions}
#'   \item{log_size}{logarithm of size of farm}
#'   \item{log_distance}{logarith of distance to nearest AHVLA center}
#'   \item{C_TYPE}{Type of farm 1 for dairy and 0 for beef}
#' }
#' @name farmsubmission
#' @references 
#' B{\"o}hning, D., Vidal‐Diez, A., Lerdsuwansri, R., Viwatwongkasem, C., and Arnold, M. (2013). "A generalization of Chao's estimator for covariate information". *Biometrics*, 69(4), 1033-1042. doi:10.1111/biom.12082
#' @usage data("farmsubmission")
"farmsubmission"


#' @title Netherlands immigrant data
#' @details Data on irregular immigrants in Netherlands
#' @docType data
#' @format Data frame with 1880 rows and 5 columns.
#' \describe{
#'   \item{capture}{Number of times a person has been captured}
#'   \item{gender}{Gender for apprehended person 1 for male and 0 for female}
#'   \item{age}{Age of apprehended person 0 for less than 40 yrs old 0 otherwise}
#'   \item{reason}{Reason for being apprehended 1 if reason was illegally
#'   residing in Netherlands 0 otherwise}
#'   \item{nation}{Nation of origin of the captured person}
#' }
#' @name netherlandsimmigrant
#' @references 
#' van Der Heijden, P. G., Bustami, R., Cruyff, M. J., Engbersen, G., and Van Houwelingen, H. C. (2003). Point and interval estimation of the population size using the truncated Poisson regression model. *Statistical Modelling*, 3(4), 305-322. doi:10.1191/1471082X03st057oa
#' @usage data("netherlandsimmigrant")
"netherlandsimmigrant"


#' @title Carcass submissions
#' @details Data on british animal farms submissions to AHVLA. British farms 
#' are able to submit samples to AHVLA if cause of death of an animal
#' cannot be determined and private veterinary surgeon decides to submit them, 
#' unless there is notifiable disease suspected then it is not required. 
#' This is data on such farms and their submissions. This
#' data set contains submission of animal carcasses only.
#' @docType data
#' @format Data frame with 1858 rows and 4 columns.
#' \describe{
#'   \item{TOTAL_SUB}{Number of submissions}
#'   \item{log_size}{logarithm of size of farm}
#'   \item{log_distance}{logarith of distance to nearest AHVLA center}
#'   \item{C_TYPE}{Type of farm 1 for dairy and 0 for beef}
#' }
#' @name carcassubmission
#' @references 
#' B{\"o}hning, D., Vidal‐Diez, A., Lerdsuwansri, R., Viwatwongkasem, C., and Arnold, M. (2013). "A generalization of Chao's estimator for covariate information". *Biometrics*, 69(4), 1033-1042. doi:10.1111/biom.12082
#' @usage data("carcassubmission")
"carcassubmission"

