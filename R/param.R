#' CC and AC parameters calculation function
#'
#' Calculates the complete-case and available-case parameters
#'
#' @param data dataset containing the variables with missing and censored values
#'
#' @return A list of length 5 containing the available-case mean and variance for each variable,
#' the complete-case covariance matrix, mean and variance for each variable.
#'
#' @export
calcu.param <- function(data) {

  if (is.null(dim(data))) stop("Need a multidimensional matrix!")

  ## CC parameters
  complete.dat <- data[complete.cases(data), ]
  CC.mean <- apply(complete.dat, 2, mean)
  CC.var <- apply(complete.dat, 2, var)
  CC.cov <- cov(complete.dat)

  ## AC parameters
  AC.mean <- numeric(ncol(data)); AC.var <- numeric(ncol(data))
  for (i in 1:ncol(data)) {
    available.dat <- na.omit(data[, i])
    AC.mean[i] <- mean(available.dat)
    AC.var[i] <- var(available.dat)
  }

  return(list(
    AC.mean = AC.mean,
    AC.var = AC.var,
    CC.cov = CC.cov,
    CC.mean = CC.mean,
    CC.var = CC.var))
}
