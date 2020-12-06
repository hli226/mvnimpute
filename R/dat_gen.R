#' Data generating function
#'
#' Generates multivariate normal data
#'
#' @param n number of observations to be generated
#' @param p number of variables to be generated
#' @param mu specified mean vector
#' @param Sigma specified covariance matrix
#' @param miss.pos a
#' @param miss.percent b
#' @param miss.type c
#' @param censor.pos d
#' @param censor.percent e
#' @param censor.type f
#'
#' @details This function generates the multivariate normal data that can be used to
#' verify the correctness of the multiple imputation algorithm. Users have to specify
#' the mean vector, variance covariance matrix, sample size and dimension of the generated
#' data, and the desired number of variables subject to missing and censoring, respectively.
#' Currently, it only supports generating the data with a certain type of MAR missing mechanism
#' and interval censoring mechanism, in which case there should be at least one variable that is
#' fully observed in the dataset.
#'
#' @export
data.generation <- function(
  n,                        # number of observations
  p,                        # number of variables
  mu,                       # mean vector
  Sigma,                    # covariance matrix
  ############################################
  miss.pos,                 # variables subject to missing
  miss.percent,             # missing percentages
  miss.type = "MCAR",       # missing mechanism
  censor.pos,               # variables subject to censoring
  censor.percent,           # censoring percentages
  censor.type = "interval"  # censoring type
) {

  #######################################
  ## Generate multivariate normal data ##
  #######################################
  ## complete data before applying the missing and censored data information
  multi.data <- mvrnorm(n = n, mu = mu, Sigma = Sigma)

  ## create a list to store the lower and upper limits of the incomplete data
  incomplete.data <- list()
  incomplete.data[[1]] <- multi.data
  incomplete.data[[2]] <- multi.data

  if (miss.type == "MCAR") {
    # gernate variable type indicator matrix
    data.ind <- MCAR.type(n, p, miss.pos, miss.percent, censor.pos, censor.percent)

    for (i in 1:n) {

      for (j in 1:p) {

        if (data.ind[i, j] == 1) {
          # observed data
          incomplete.data[[1]][i, j] <- incomplete.data[[2]][i, j] <- multi.data[i, j]
        } else if (data.ind[i, j] == 0) {
          # missing data
          incomplete.data[[1]][i, j] <- -10e10 # a small number for negative infinity
          incomplete.data[[2]][i, j] <- 10e10  # a large number for positive infinity
        } else {
          # censored data
          if (censor.type == "interval") {
            # interval censoring
            c1 <- rnorm(1, 0, 1)
            c2 <- rnorm(1, 2, 1)

            incomplete.data[[1]][i, j] <- min(c1, c2)
            incomplete.data[[2]][i, j] <- max(c1, c2)
          } else if (censor.type == "right") {
            # right censoring
            c1 <- rnorm(1, 0, 1)

            incomplete.data[[1]][i, j] <- c1
            incomplete.data[[2]][i, j] <- -10e10 # a small number for negative infinity
          } else if (censor.type == "left"){
            # left censoring
            c2 <- rnorm(1, 2, 1)

            incomplete.data[[1]][i, j] <- c2
            incomplete.data[[2]][i, j] <- 10e10 # a large number for positive infinity
          }
        }
      }
    }
    # rename each variable
    colnames(incomplete.data[[1]]) <- colnames(incomplete.data[[2]]) <- colnames(data.ind)
  }
  else if (miss.type == "MAR") {

    # observed data and missing data
    incomplete.data <- MAR.type(incomplete.data, miss.pos, miss.percent, censor.pos, censor.percent)$data
    data.ind <- MAR.type(incomplete.data, miss.pos, miss.percent, censor.pos, censor.percent)$ind

    # censored data
    if (censor.type == "interval") {
      for (i in 1:length(censor.pos)) {

        # the observations that are not missing
        not.miss <- which(incomplete.data[[1]][, censor.pos[i]] > -10e10)

        # introduce randomness into the censoring percentage
        c.percent <- censor.percent[i] * runif(1, 0.95, 1.05)
        censor.obs <- numeric(n * c.percent)

        for (j in 1:length(censor.obs)) {
          censor.obs <- sample(not.miss, n * c.percent, replace = FALSE)
        }
        c1 <- rnorm(n * c.percent, 0, 1)
        c2 <- rnorm(n * c.percent, 2, 1)

        data.ind[censor.obs, censor.pos[i]] <- 2
        incomplete.data[[1]][censor.obs, censor.pos[i]] <- pmin(c1, c2)
        incomplete.data[[2]][censor.obs, censor.pos[i]] <- pmax(c1, c2)
      }
      data.ind <- ifelse(is.na(data.ind), 1, data.ind)
    }
    else if (censor.type == "right") {
      for (i in 1:length(censor.pos)) {

        # the observations that are not missing
        not.miss <- which(incomplete.data[[1]][, censor.pos[i]] > -10e10)

        # introduce randomness into the censoring percentage
        c.percent <- censor.percent[i] * runif(1, 0.95, 1.05)
        censor.obs <- numeric(n * c.percent)

        for (j in 1:length(censor.obs)) {
          censor.obs <- sample(not.miss, n * c.percent, replace = FALSE)
        }
        c1 <- rnorm(n * c.percent, 0, 1)

        data.ind[censor.obs, censor.pos[i]] <- 2
        incomplete.data[[1]][censor.obs, censor.pos[i]] <- pmin(c1, c2)
        incomplete.data[[2]][censor.obs, censor.pos[i]] <-10e10
      }
      data.ind <- ifelse(is.na(data.ind), 1, data.ind)
    } else if (censor.type == "left") {
      for (i in 1:length(censor.pos)) {

        # the observations that are not missing
        not.miss <- which(incomplete.data[[1]][, censor.pos[i]] > -10e10)

        # introduce randomness into the censoring percentage
        c.percent <- censor.percent[i] * runif(1, 0.95, 1.05)
        censor.obs <- numeric(n * c.percent)

        for (j in 1:length(censor.obs)) {
          censor.obs <- sample(not.miss, n * c.percent, replace = FALSE)
        }
        c2 <- rnorm(n * c.percent, 2, 1)

        data.ind[censor.obs, censor.pos[i]] <- 2
        incomplete.data[[1]][censor.obs, censor.pos[i]] <- -10e10
        incomplete.data[[2]][censor.obs, censor.pos[i]] <- pmax(c1, c2)
      }
      data.ind <- ifelse(is.na(data.ind), 1, data.ind)
    }
  }

  return(list(data = incomplete.data, ind = data.ind))

}
