#' Initial data fill-in fucntion
#'
#' Fill in data with missing and censored values for multiple imputation
#'
#' @param dat dataset with missing and censored values
#' @param method method to be used for initial filled-in values
#'
#' @return a complete dataset with missing values filled in
#' @export
fill.dat <- function(dat, method = "normal")
{
  fill.dat <- matrix(nrow = nrow(dat), ncol = ncol(dat))
  if (method == "normal") {
    for (i in 1:ncol(dat)) {
      m <- mean(na.omit(dat[, i])); v <- var(na.omit(dat[, i]))
      fill.dat[, i] <- ifelse(is.na(dat[, i]), rnorm(1, m, sqrt(v)), dat[, i])
    }
  }
  return(fill.dat)
}
