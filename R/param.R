#' CC and AC parameters calculation function
#'
#' Calculates the complete-case and avaliable-case parameters
#'
#' @param dat dataset containing the variables with missing and censored values.
#'
#' @return a list of length 5 containing the avaliable-case mean and variance for each variable, the complete-case covariance matrix,
#' the complete-case mean and variance for each variable.
#'
#' @export
calcu.param <- function(dat) {

  if (is.null(dim(dat))) stop("Need a multidimensional matrix!")

  ## CC parameters
  complete.dat <- dat[complete.cases(dat), ]
  CC.mean <- apply(complete.dat, 2, mean)
  CC.var <- apply(complete.dat, 2, var)
  CC.cov <- cov(complete.dat)

  ## AC parameters
  AC.mean <- numeric(ncol(dat)); AC.var <- numeric(ncol(dat))
  for (i in 1:ncol(dat)) {
    available.dat <- na.omit(dat[, i])
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
