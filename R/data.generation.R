#' Data generating function
#'
#' Generates multivariate normal data with missing and censored values
#'
#' @param n number of observations to be generated
#' @param p number of variables to be generated
#' @param mu specified mean vector
#' @param Sigma specified covariance matrix
#' @param miss.pos indexes of the variables that have missing values
#' @param miss.percent the missing percentage in each variable that has missing values
#' @param miss.type mechanism of missing data, available options include the default \code{"MCAR"} for missing completely at random, and
#' \code{"MAR"} for missing at random
#' @param censor.pos indexes of the variables that have censored values
#' @param censor.percent1 the percentage of the lower limit of the censored values
#' @param censor.percent2 the percentage of the upper limit of the censored values
#' @param censor.type type of censored data, available options include the default \code{"interval"} for interval censoring, \code{"right"}
#' for right censoring, and \code{"left"} for left censoring
#'
#' @details This function generates the multivariate normal data  with missing and censored values, that can be used to verify the
#' correctness of the multiple imputation algorithm. Users have to specify the sample size \code{n}, the dimension \code{p}, the
#' mean vector \code{mu}, and the covariance matrix \code{Sigma} of the data. Additionally, the information for the missing and censoring,
#' such as the indexes, percentages for missing and censored values, missing mechanisms or censoring types, should also be included.
#' Specifically, we assign a pair of values, say \code{X_{ll}} and \code{X_{ul}}, to each of the data point in the following way: if the
#' values is observed, we set \code{X_{ll}} = \code{X_{ul}}; if the value is missing, we set \code{X_{ll}} = -10e10 and \code{X_{ul}} =
#' 10e10; if the value is censored, we set \code{X_{ll}} and \code{X_{ll}} to finite values according to the specified censoring type.
#'
#' @return A list of two containing two matrices.
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
  censor.percent1,          # censoring percentages
  censor.percent2,
  censor.type = "interval"  # censoring type
) {

  ## generate data
  dat <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  data.ind <- matrix(NA, nrow = nrow(dat), ncol = ncol(dat))

  incomplete.data <- list()
  incomplete.data[[1]] <- incomplete.data[[2]] <- matrix(NA, nrow = nrow(dat), ncol = ncol(dat))
  ##### missing values
  ## MCAR
  if (miss.type == "MCAR" ) {

    for (i in 1:length(miss.pos)) {
      m.percent <- miss.percent[i]
      data.ind[, miss.pos[i]] <- sample(c(0, 1), n, replace = TRUE, prob = c(m.percent, 1 - m.percent))
    }

    # set all other entries to 1
    data.ind <- ifelse(is.na(data.ind), 1, data.ind)

    # dataset
    incomplete.data[[1]] <- ifelse(data.ind == 0, -10e10, dat)
    incomplete.data[[2]] <- ifelse(data.ind == 0, 10e10, dat)

    ##### censored values

    for (i in 1:length(censor.pos)) {

      not.miss <- which(data.ind[, censor.pos[i]] != 0)
      data.remain <- dat[not.miss, censor.pos[i]]
      data.remain.sort <- sort(data.remain)

      ## interval censoring
      if (censor.type == "interval") {
        # limits
        ll <- data.remain.sort[round(length(data.remain) * censor.percent1[i])]
        ul <- data.remain.sort[round(length(data.remain) * censor.percent2[i])]
        for (j in 1:length(data.remain)) {
          if (data.remain[j] >= ll & data.remain[j] <= ul) {
            incomplete.data[[1]][not.miss, censor.pos[i]][j] <- ll
            incomplete.data[[2]][not.miss, censor.pos[i]][j] <- ul
            # updat data type indicator matrix
            data.ind[not.miss, censor.pos[i]][j] <- 2
          }
        }
      }
      ## right censoring
      else if (censor.type == "right") {
        if (censor.percent2[i] != 10e10) stop("censor.percent2 should be infinity in right censoring!")
        # limit
        ll <- data.remain.sort[round(length(data.remain) * censor.percent1[i])]
        for (j in 1:length(data.remain)) {
          if (data.remain[j] >= ll) {
            incomplete.data[[1]][not.miss, censor.pos[i]][j] <- ll
            incomplete.data[[2]][not.miss, censor.pos[i]][j] <- 10e10
            # updat data type indicator matrix
            data.ind[not.miss, censor.pos[i]][j] <- 2
          }
        }
      }
      ## left censoring
      else if (censor.type == "left") {
        if (censor.percent1[i] != -10e10) stop("censor.percent2 should be negative infinity in left censoring!")
        # limit
        ul <- data.remain.sort[round(length(data.remain) * censor.percent2[i])]
        for (j in 1:length(data.remain)) {
          if (data.remain[j] <= ul) {
            incomplete.data[[1]][not.miss, censor.pos[i]][j] <- -10e10
            incomplete.data[[2]][not.miss, censor.pos[i]][j] <- ul
            # updat data type indicator matrix
            data.ind[not.miss, censor.pos[i]][j] <- 2
          }
        }
      }
    }
  }
  ## MAR
  else if (miss.type == "MAR") {
    # the index of the observed variable that the missing variable depend on
    total.var <- 1:p
    not.miss <- total.var[!total.var %in% miss.pos]
    obs.pos <- not.miss[!not.miss %in% censor.pos][1]

    for (i in 1:length(miss.pos)) {

      ## missing percentage
      m.percent <- miss.percent[i]
      cutoff <- sort(dat[, miss.pos[i]])[round(n * m.percent)]

      data.ind[, miss.pos[i]] <- ifelse(dat[, obs.pos] <= cutoff, 0, 1)
    }
    # set all other entries to 1
    data.ind <- ifelse(is.na(data.ind), 1, data.ind)

    # dataset
    incomplete.data[[1]] <- ifelse(data.ind == 0, -10e10, dat)
    incomplete.data[[2]] <- ifelse(data.ind == 0, 10e10, dat)

    ##### censored values
    for (i in 1:length(censor.pos)) {

      not.miss <- which(data.ind[, censor.pos[i]] != 0)
      data.remain <- dat[not.miss, censor.pos[i]]
      data.remain.sort <- sort(data.remain)

      ## interval censoring
      if (censor.type == "interval") {
        # limits
        ll <- data.remain.sort[round(length(data.remain) * censor.percent1[i])]
        ul <- data.remain.sort[round(length(data.remain) * censor.percent2[i])]
        for (j in 1:length(data.remain)) {
          if (data.remain[j] >= ll & data.remain[j] <= ul) {
            incomplete.data[[1]][not.miss, censor.pos[i]][j] <- ll
            incomplete.data[[2]][not.miss, censor.pos[i]][j] <- ul
            # updat data type indicator matrix
            data.ind[not.miss, censor.pos[i]][j] <- 2
          }
        }
      }
      ## right censoring
      else if (censor.type == "right") {
        if (censor.percent2[i] != 10e10) stop("censor.percent2 should be infinity in right censoring!")
        # limit
        ll <- data.remain.sort[round(length(data.remain) * censor.percent1[i])]
        for (j in 1:length(data.remain)) {
          if (data.remain[j] >= ll) {
            incomplete.data[[1]][not.miss, censor.pos[i]][j] <- ll
            incomplete.data[[2]][not.miss, censor.pos[i]][j] <- 10e10
            # update data type indicator matrix
            data.ind[not.miss, censor.pos[i]][j] <- 2
          }
        }
      }
      ## left censoring
      else if (censor.type == "left") {
        if (censor.percent1[i] != -10e10) stop("censor.percent1 should be negative infinity in right censoring!")
        # limit
        ul <- data.remain.sort[round(length(data.remain) * censor.percent2[i])]
        for (j in 1:length(data.remain)) {
          if (data.remain[j] <= ul) {
            incomplete.data[[1]][not.miss, censor.pos[i]][j] <- -10e10
            incomplete.data[[2]][not.miss, censor.pos[i]][j] <- ul
            # updat data type indicator matrix
            data.ind[not.miss, censor.pos[i]][j] <- 2
          }
        }
      }
    }
  }
  colnames(incomplete.data[[1]]) <- colnames(incomplete.data[[2]]) <- colnames(data.ind) <- colnames(dat) <- paste0("y", 1:p)
  return(list(full.data = dat, data = incomplete.data, indicator = data.ind))

}
