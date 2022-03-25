#' @title All farm submissions
#' @details Data on british animal farms submissions to AHVLA. It is required
#' for british farms to submit samples to AHVLA if unknown disease
#' if private veterinary surgeons cannot determine the cause of
#' death for the animal and decide to submit them, unless there is notifiable
#' disease suspected then it is required. This is data on such farms. This
#' data set contains farms that have submitted anything at all not only carcasses but
#' also blood samples etc.
#' @docType data
#' @format Data frame with 12036 rows and 4 columns.
#' \describe{
#'   \item{TOTAL_SUB}{Number of submissions}
#'   \item{log_size}{logarithm of size of farm}
#'   \item{log_distance}{logarith of distance to nearest AHVLA center}
#'   \item{C_TYPE}{Type of farm 1 for dairy and 0 for beef}
#' }
#' @name farmsubmission
#' @references Bohning et.al (2013) BIOMETRICS 69, 1033-1042
#' @usage data(farmsubmission)
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
#' @references van der Heijden et.al (2003) Statistical Modelling 3, 305-322
#' @usage data(netherlandsimmigrant)
"netherlandsimmigrant"
#' @title Carcass submissions
#' @details Data on british animal farms submissions to AHVLA. It is required
#' for british farms to submit samples to AHVLA if unknown disease
#' if private veterinary surgeons cannot determine the cause of
#' death for the animal and decide to submit them, unless there is notifiable
#' disease suspected then it is required. This is data on such farms. This
#' data set contains farms that have submitted animal carcasses.
#' @docType data
#' @format Data frame with 1858 rows and 4 columns.
#' \describe{
#'   \item{TOTAL_SUB}{Number of submissions}
#'   \item{log_size}{logarithm of size of farm}
#'   \item{log_distance}{logarith of distance to nearest AHVLA center}
#'   \item{C_TYPE}{Type of farm 1 for dairy and 0 for beef}
#' }
#' @name carcassubmission
#' @references Bohning et.al (2013) BIOMETRICS 69, 1033-1042
#' @usage data(carcassubmission)
"carcassubmission"
