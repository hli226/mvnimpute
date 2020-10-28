#' Data generating function
#'
#' Generates multivariate normal data.
#'
#' @param n number of observations to be generated.
#' @param p number of variables to be generated.
#' @param m.vec specified mean vector.
#' @param sig specified covariance matrix.
#' @param miss.num number of variables having missing values.
#' @param censor.num number of variables having censored values.
#'
#' @details The function generates the multivariate normal data that can be used to verify the correctness of the multiple imputation
#' algorithm. Users have to specify the mean vector, variance covariance matrix, sample size and dimension of the generated data, and
#' the desired number of variables subject to missing and censoring, respectively. Currently, it only supports generating the data with
#' a certain type of MAR missing mechanism and interval censoring mechanism, in which case there should be at least one variable that is fully
#' observed in the dataset.
#'
#' @examples
#' n <- 1000
#' p <- 3
#' m.vec <- c(3, 2, 1)
#' V <- diag(p)
#' miss.num <- 1
#' censor.num <- 1
#' dat <- dat.gen(n, m.vec, V, p, miss.num, censor.num)
#'
#' @return A list containing the full data, the complete data with missing and censoring information applied,
#'   the missing and censoring index, positions of the variables including the missing or censored values, or
#'   the fully observed variable, and the two limits of the censored values
#'
#' @export
dat.gen <- function(
  n,          # number of data to be generated
  p,          # number of in the variables of the data
  m.vec,      # mean vector should have length p
  sig,        # variance-covariance matrix have dimension p * p
  miss.num,   # number of variables subject to missing
  censor.num  # number of censoring subject
){

  `%notin%` <- Negate(`%in%`)

  if (length(m.vec) != p) stop("Your desired dimensions do not match!")
  if (dim(sig)[1] != dim(sig)[2]) stop("Covariance matrix should be a sqaure matrix!")

  # use rmvnorm function from the mvtnorm package to generate multivariate normal data
  dat <- rmvnorm(n, m.vec, sig)

  # full data without applying missing and censoring information
  full <- dat

  ### generate missing and censoring index
  # will draw random number for the missing variables and censoring variables
  var.num <- 1:p

  # 1. missing data and censoring data in the dataset
  # (1) when data has a mixed type of missing and censoring
  if (miss.num != 0 & censor.num != 0)
  {
    miss.pos <- sample(var.num, miss.num, replace = FALSE)
    censor.pos <- sample(var.num[var.num %notin% miss.pos], censor.num, replace = FALSE)
    full.pos <- var.num[var.num %notin% c(miss.pos, censor.pos)]

    # fully observed data that missing and censoring data will depend on
    # the first of the fully observed variables
    fully.obs <- as.matrix(dat[, full.pos], nrow = n, byrow = FALSE)[, 1]

    ### missing data
    # generate missing index
    miss.indx <- matrix(NA, nrow = n, ncol = miss.num)

    miss.point <- numeric(miss.num)

    for (i in 1:miss.num) {

      # randomly generate truncation point
      miss.prob <- sample(seq(0.1, 0.3, by = 0.001), 1, replace = TRUE)
      miss.point[i] <- sort(fully.obs)[miss.prob * length(fully.obs)]

      # MAR missing mechanism
      miss.indx[, i]  <- ifelse(fully.obs < miss.point[i], 1, 0)

      # set the value where missingness happens to be NA
      miss.p <- miss.pos[i];  miss.in <- miss.indx[, i]
      dat[, miss.p][miss.in == 1] <- NA
    }

    ### censoring data
    ### generate censoring times

    censor.ll <- matrix(NA, nrow = n, ncol = censor.num) # lower limit of the censoring interval
    censor.ul <- matrix(NA, nrow = n, ncol = censor.num) # upper limit of the censoring interval
    censor.indx <- matrix(NA, nrow = n, ncol = censor.num) # censoring indices
    up <- quantile(fully.obs, prob = 0.45)
    down <- quantile(fully.obs, prob = 0.65)

    for (i in 1:censor.num) {

      censor.p <- censor.pos[i]

      # censoring times

      t1 <- rexp(n, 1 / up)
      t2 <- rexp(n , 1 / down)
      censor.ll[, i] <- pmin(t1, t2)
      censor.ul[, i] <- pmax(t1, t2)

      # generating censoring indices

      censor.dat <- dat[, censor.p]
      censor.indx[, i] <- ifelse((censor.dat > censor.ll[, i] & censor.dat < censor.ul[, i]), 1, 0)

      # set the value where data is censoring to NaN

      censor.p <- censor.pos[i]; censor.in <- censor.indx[, i]
      dat[, censor.p][censor.in == 1] <- NaN
    }

    # rename the columns
    # the prefix suggests the type of the variable, and the number indicates the column number
    colnames(dat) <- 1:p
    colnames(dat)[miss.pos] <- paste("miss", miss.pos, sep = ".")
    colnames(dat)[censor.pos] <- paste("censor", censor.pos, sep = ".")

    f <- 1:p
    m.s <- c(miss.pos, censor.pos)
    full.f <- f[-m.s]
    colnames(dat)[-c(miss.pos, censor.pos)] <- paste("full", full.f, sep = ".")

    # create a list to store the results
    res <- list(
      full.dat = full,
      comp.dat = dat,
      miss.indx = miss.indx, miss.pos = miss.pos, miss.point = miss.point,
      censor.indx = censor.indx, censor.pos = censor.pos,
      C1 = censor.ll, C2 = censor.ul,
      full.pos = full.pos
    )
  }

  # (2) only missing data in the dataset

  if (miss.num != 0 & censor.num == 0)
  {
    miss.pos <- sample(var.num, miss.num, replace = FALSE)
    full.pos <- var.num[var.num %notin% miss.pos]

    # fully observed data that missing data will depend on
    # the first of the fully observed variables
    fully.obs <- as.matrix(dat[, full.pos], nrow = n, byrow = FALSE)[, 1]

    ### missing data
    # generate missing index
    miss.indx <- matrix(NA, nrow = n, ncol = miss.num)
    miss.point <- numeric(miss.num)

    for (i in 1:miss.num) {

      # randomly generate missing probability
      miss.prob <- sample(seq(0.1, 0.3, by = 0.001), 1, replace = TRUE)
      miss.point[i] <- sort(fully.obs)[miss.prob * length(fully.obs)]

      # MAR missing mechanism
      miss.indx[, i]  <- ifelse(fully.obs < miss.point[i], 1, 0)

      miss.p <- miss.pos[i];  miss.in <- miss.indx[, i]

      # set the value where data is missing to NA
      dat[, miss.p][miss.in == 1] <- NA
    }

    # rename the columns
    colnames(dat) <- 1:p
    colnames(dat)[miss.pos] <- paste("miss", miss.pos, sep = ".")

    f <- 1:p
    m.s <- miss.pos
    full.f <- f[-m.s]
    colnames(dat)[-miss.pos] <- paste("full", full.f, sep = ".")

    # create a list to store the results
    res <- list(
      full.dat = full,
      comp.dat = dat,
      miss.indx = miss.indx, miss.pos = miss.pos, miss.point = miss.point,
      full.pos = full.pos
    )
  }

  # (3) only censoring data in the dataset

  if (miss.num == 0 & censor.num != 0)
  {
    censor.pos <- sample(var.num, censor.num, replace = FALSE)
    full.pos <- var.num[var.num %notin% censor.pos]
    fully.obs <- as.matrix(dat[, full.pos], nrow = n, byrow = FALSE)[, 1] # fully observed data that censoring data will depen

    ### censoring data
    ### generate censoring times
    censor.ll <- matrix(NA, nrow = n, ncol = censor.num)
    censor.ul <- matrix(NA, nrow = n, ncol = censor.num)
    censor.indx <- matrix(NA, nrow = n, ncol = censor.num)
    up <- quantile(fully.obs, prob = 0.45)
    down <- quantile(fully.obs, prob = 0.65)

    for (i in 1:censor.num) {
      censor.p <- censor.pos[i]
      t1 <- rexp(n, 1 / up)
      t2 <- rexp(n , 1 / down)

      censor.ll[, i] <- pmin(t1, t2)
      censor.ul[, i] <- pmax(t1, t2)

      censor.dat <- dat[, censor.p]

      censor.indx[, i] <- ifelse((censor.dat > censor.ll[, i] & censor.dat < censor.ul[, i]), 1, 0)

      censor.p <- censor.pos[i]; censor.in <- censor.indx[, i]
      dat[, censor.p][censor.in == 1] <- NaN
    }

    # rename the columns
    colnames(dat) <- 1:p
    colnames(dat)[censor.pos] <- paste("censor", censor.pos, sep = ".")

    f <- 1:p
    m.s <- censor.pos
    full.f <- f[-m.s]
    colnames(dat)[-censor.pos] <- paste("full", full.f, sep = ".")

    # create a list to store the results
    res <- list(
      full.dat = full,
      comp.dat = dat,
      censor.indx = censor.indx, censor.pos = censor.pos,
      C1 = censor.ll, C2 = censor.ul,
      full.pos = full.pos
    )
  }

  return(res)
}
