#' Function to calculate complete case and avaliable case parameters
#'
#' @param dat dataset to compute the complete case and available case parameters
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
