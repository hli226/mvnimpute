#' Combined NHANES dataset from 2011 - 2016
#'
#' A dataset containing the demographic, PCB and metal pollutant measurements variables.
#'
#' @format A data set frame with 14268 rows and 154 variables:
#' \describe{
#'   \item{SEQN}{respondent sequence nubmer}
#'   \item{SDDSRVYR}{survey cycle identifier, it takes 7, 8, 9 in the data}
#'   \item{WTINT2YR}{full data 2 year interview weight}
#'   \item{WTMEC2YR}{full data 2 year MEC exam weight}
#'   \item{SDMVPSU}{Masked variance unit pseudo-PSU variable for variance estimation}
#'   \item{SDMVSTRA}{Masked variance unit pseudo-stratum variable for variance estimation}
#' }
#'
#' @details The dataset is combined from the NHANES release cycle 2011-2012, 2013-2014, 2015-2016. It includes three domains, the
#' demographics data, PCB and metal pollutant lab measurements data. There are lots of missing and censored values in this dataset.
#' A value is left censored if it falls below the corresponding lower detection of limit (LLOD). Originally, the censored values were
#' imputed as the LLOD divided by the square root of 2. For PCB lab measurements data, a pooled sample design was used to minimize
#' the potential harms it may have to the survey participants. Participants are divided into several pools based on their demographic
#' status, and within each sample pool, the pooled lab measurements were replicated for each individual.
#'
#' @source \url{https://www.cdc.gov/nchs/nhanes/index.htm}
"nhanes.dat"
