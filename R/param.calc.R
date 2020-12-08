#' CC and AC parameters calculation function
#'
#' Calculates the complete-case and available-case parameters
#'
#' @param data data set containing the missing and censored values
#'
#' @return A list of four containing the available-case mean for each variable,
#' the complete-case covariance matrix, the available-case mean and variance for each variable.
#'
#' @export
param.calc <- function(data) {

  lvalue <- data[[1]]
  rvalue <- data[[2]]

  # number of observations and variables
  n <- nrow(lvalue); p <- ncol(rvalue)

  ## complete-case data
  CC.dat <- lvalue[apply(lvalue == rvalue, 1, sum) == p, ]

  CC.mean <- apply(CC.dat, 2, mean)
  CC.cov <- cov(CC.dat)

  ## available-case dataa: does not include censored values
  AC.mean <- numeric(p)
  AC.var <- numeric(p)
  names(AC.mean) <- names(AC.var) <- colnames(lvalue)
  for (i in 1:p) {

    obs.indx <- lvalue[, i] == rvalue[, i]
    obs.dat <- lvalue[obs.indx, i]

    AC.mean[i] <- mean(obs.dat)
    AC.var[i] <- var(obs.dat)
  }

  return(list(
    CC.mean = CC.mean,   # complete-case mean vector
    CC.cov = CC.cov,     # complete-case covariance matrix
    AC.mean = AC.mean,   # available-case mean vector
    AC.var = AC.var      # available-case variance
  ))
}
