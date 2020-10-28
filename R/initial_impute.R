#' Initial missing data fill-in function
#'
#' Fills in data with missing and censored values for multiple imputation
#'
#' @param data dataset with missing and censored values
#' @param miss.index matrix containing the missing index
#' @param miss.pos position of variables with missing values in the original dataset
#' @param censor.index matrix containing the censoring index
#' @param censor.pos position of variables with censored values in the original dataset
#'
#' @details The unobserved values will initially be filled in by some single imputation methods.
#' Currently, it only supports generating random values from the normal distribution with the
#' mean and variance as the complete-case mean and variance, respectively.
#'
#' @return a complete dataset with missing values filled in
#' @export
initial.impute <- function(data,
                           miss.index = NULL, miss.pos = NULL,
                           censor.index = NULL, censor.pos = NULL
)
{
  miss.dat <- data

  # fill-in missing data

  if (!is.null(miss.index)){
    for (i in 1:length(miss.pos)) {

      miss <- miss.pos[i]
      if (!is.null(dim(miss.index))) miss.in <- miss.index[, i]
      else miss.in <- miss.index

      x <- miss.dat[, miss][miss.in == 1]
      miss.mean <- calcu.param(miss.dat)$CC.mean[miss]
      miss.var <- calcu.param(miss.dat)$CC.var[miss]

      for (j in 1:length(x)) {
        x[j] <- rnorm(1, miss.mean, sqrt(miss.var))
        # generate normal random variable using the complete cases mean and variance
      }

      message(paste("Impute missing data column ", i, sep = ""))

      miss.dat[, miss][miss.in == 1] <- x

    }
  }


  # fill-in censored data
  if (!is.null(censor.index)) {
    for (i in 1:length(censor.pos)) {

      censor <- censor.pos[i]
      if (!is.null(dim(censor.index))) censor.in <- censor.index[, i]
      else censor.in <- censor.index

      t. <- miss.dat[, censor][censor.in == 1]
      censor.mean <- calcu.param(miss.dat)$CC.mean[censor]
      censor.var <- calcu.param(miss.dat)$CC.var[censor]

      for (j in 1:length(t.)) {

        t.[j] <- rnorm(1,
                       mean = censor.mean, sqrt(censor.var))
        # generate random variable uniformly in the interval
      }

      message(paste("Impute censored data column ", i, sep = ""))

      miss.dat[, censor][censor.in == 1] <- t.

    }
  }

  return(as.matrix(miss.dat))
}
