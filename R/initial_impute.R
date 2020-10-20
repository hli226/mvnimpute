#' Initial missing data fill-in fucntion
#'
#' Fills in data with missing and censored values for multiple imputation
#'
#' @param dat dataset with missing and censored values
#' @param miss.indx matrix containing the missing index
#' @param miss.pos position of variables with missing values in the original dataset
#' @param censor.indx matrix containing the censoring index
#' @param censor.pos position of variables with censored values in the original dataset
#'
#' @details The missing values will initially be filled in by some single imputation methods. Currently, it only supports generating
#' random values from the normal distribution with the mean and variance as the complete-case mean and variance, respectively.
#'
#' @return a complete dataset with missing values filled in
#' @export
initial.impute <- function(dat, # input dataset, it should be data with missing info
                           miss.indx, miss.pos,
                           censor.indx, censor.pos
)
{
  miss.dat <- dat

  # dat <- miss.dat[complete.cases(miss.dat), ]

  # fill-in missing data

  if (!is.null(miss.indx)){
    for (i in 1:length(miss.pos)) {

      miss <- miss.pos[i]; miss.in <- miss.indx[, i]
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
  if (!is.null(censor.indx)) {
    for (i in 1:length(censor.pos)) {

      censor <- censor.pos[i]; censor.in <- censor.indx[, i]
      t. <- miss.dat[, censor][censor.in == 1]
      censor.mean <- calcu.param(miss.dat)$CC.mean[censor]
      censor.var <- calcu.param(miss.dat)$CC.var[censor]
      # censor.l <- censor.val[[1]][censor.in == 1]; censor.u <- censor.val[[2]][censor.in == 1]

      for (j in 1:length(t.)) {

        # t.[j] <- rtruncnorm(1, a = censor.l[j], b = censor.u[j],
        #                     mean = censor.mean, sqrt(censor.var))

        # t.[j] <- rnorm(1,
        #                mean = censor.mean, sqrt(censor.var))
        t.[j] <- rnorm(1,
                       mean = censor.mean, sqrt(censor.var))
        # generate random variable uniformly in the interval
      }

      message(paste("Impute censored data column ", i, sep = ""))

      miss.dat[, censor][censor.in == 1] <- t.

    }
  }

  return(miss.dat)
}
