#' Initial missing data fill-in function
#'
#' Fills in data with missing and censored values for multiple imputation
#'
#' @param data dataset with missing and censored values
#'
#' @details The unobserved values are initially filled in by some single imputation methods.
#' Currently, it only supports generating random values from the normal distribution with the
#' mean and variance as the complete-case mean and variance, respectively.
#'
#' @return a complete dataset with missing and censored values filled in
#' @export
#### simple imputation function to make up incomplete data
single.imputation <- function(
  data      # the input dataset
) {

  CC.mean <- param.calc(data)$CC.mean
  CC.var <- diag(param.calc(data)$CC.cov)

  lvalue <- data[[1]]
  rvalue <- data[[2]]

  # number of observations and variables
  n <- nrow(lvalue); p <- ncol(lvalue)

  # create a matrix for storing the resulting data from simple imputation that makes up the incomplete data
  sim.imp.data <- matrix(NA, nrow = n, ncol = p)

  for (i in 1:n) {

    for (j in 1:p) {

      # (1) observed values
      if (lvalue[i, j] == rvalue[i, j]) {
        sim.imp.data[i, j] <- lvalue[i, j]
      }

      # (2) missing values
      else if (lvalue[i, j] == -10e10 & rvalue[i, j] == 10e10) {
        sim.imp.data[i, j] <- rnorm(1, mean = CC.mean[j], sd = sqrt(CC.var[j]))
      }

      # (3) censored values
      else if (lvalue[i, j] > -10e10 & rvalue[i, j] < 10e10) { # interval censoring
        sim.imp.data[i, j] <- rtruncnorm(1, a = lvalue[i, j], b = rvalue[i, j],
                                         mean = CC.mean[j], sd = sqrt(CC.var[j]))
      }
      else if (lvalue[i, j] > -10e10 & rvalue[i, j] == 10e10) { # right censoring
        sim.imp.data[i, j] <- rtruncnorm(1, a = lvalue[i, j],
                                         mean = CC.mean[j], sd = sqrt(CC.var[j]))
      }
      else if (lvalue[i, j] == -10e10 & rvalue[i, j] < 10e10) { # left censoring
        sim.imp.data[i, j] <- rtruncnorm(1, b = rvalue[i, j],
                                         mean = CC.mean[j], sd = sqrt(CC.var[j]))
      }
    }
  }

  # rename the columns of the imputed data
  colnames(sim.imp.data) <- colnames(lvalue)

  return(sim.imp.data)
}
