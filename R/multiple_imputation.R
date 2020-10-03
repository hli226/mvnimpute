#' Multiple imputation function
#'
#' Implements the multiple imputation of the missing and censored values
#'
#' @param iter number of iterations
#' @param prior.param list of prior parameter values
#' @param ini.vals list of initial values
#' @param dat dataset to do multiple imputation
#' @param miss.indx matrix of missing index
#' @param miss.pos vector containing the positions of variables with missing values in the original dataset
#' @param censor.indx matrix of censoring index
#' @param censor.pos vector containing the position of variables with censored values in the original dataset
#' @param censor.val list containing cutoff values for the censored variables
#' @param censor.type discrete variable specifying the censoring type, it takes "left" for left-censoring, "right" for right-censoring,
#' and "interval" for interval-censoring
#'
#' @details This function implements the multiple imputation algorithm that concurrently handles missing and censored values. This
#' function requires matrices that contain the respective indices for missing or censored values, and vectors that specify the
#' positions of the variables having the missing or censored values in the original dataset. \code{Censor.val} includes the cutoff
#' values for the censored values, it should be a matrix which has the same dimension as the censoring index matrix if the data
#' is right- or left-censored, while it should be a list of length 2 if the data is interval-censored. The first element of the
#' list should contain the lower bounds of the censored values, and the upper bounds should be in the second element. \code{prior.param}
#' is a list containing the specified parameters of the prior distributions, and \code{ini.vals} is a list containing the chosen
#' inital values for the parameters.
#'
#' @return A list of length 4 containing the simulated mean and variance vectors, covariance matrix and the imputed dataset from
#' each iteration.
#'
#' @references
#' Tanner, M., & Wong, W. (1987). The Calculation of Posterior Distributions by Data Augmentation.
#' \emph{Journal of the American Statistical Association}, \bold{82(398)}, 528-540.
#'
#' @export
multiple.impute <- function(
  iter,                     # number of iterations to run
  prior.param,              # list of prior parameters
  ini.vals,                 # list of initial values
  dat,                      # dataset to do imputation (with missing values)
  miss.indx,                # matrix of  missing index
  miss.pos,
  censor.indx,              # matrix of censoring index
  censor.pos,
  censor.val,
  censor.type = "interval"
)
{

  iter.dat <- dat

  # fill.dat. <- initial.impute(dat, miss.pos, censor.pos)

  n <- nrow(iter.dat); p <- ncol(iter.dat)
  # fill the missing values

  # prior parameters
  ### prior specification
  mu.0 <- prior.param[[1]]
  Lambda.0 <- prior.param[[2]]
  kappa.0 <- prior.param[[3]]
  nu.0 <- prior.param[[4]]

  # initial values
  mu.iter <- mu.ini <- ini.vals[[1]]
  sig.iter <- sig.ini <- ini.vals[[2]]

  # vector and list to store results
  impute <- list()
  Mu.iter <- matrix(nrow = iter + 1, ncol = ncol(iter.dat))
  Sig.iter <- matrix(nrow = iter + 1, ncol = ncol(iter.dat))
  Covmat <- list()
  Mu.iter[1, ] <- mu.ini
  Sig.iter[1, ] <- diag(sig.ini)
  Covmat[[1]] <- sig.ini
  cond. <- list()

  colnames(Mu.iter) <- colnames(Sig.iter) <- colnames(dat)

  # posterior parameters that do not depend on the data
  kappa.n <- kappa.0 + n
  nu.n <- nu.0 + n

  for (i in 1:iter) { ## iteration number

    ### posterior parameters
    y.bar <- apply(iter.dat, 2, mean)
    mu.n <- kappa.0 * mu.0 / (kappa.0 + n) + n  * y.bar / (kappa.0 + n)

    # P-step
    # update mu vector from normal distribution condition on Sigma
    mu.iter <- rmvnorm(1, mean = mu.n, sig.iter/kappa.n)

    # posterior parameters for covariance matrix
    S <- apply(iter.dat, 1, "-", mu.iter) %*%
      t(apply(iter.dat, 1, "-", mu.iter))

    Lambda.n <- Lambda.0 + S + kappa.0 * n * (y.bar - mu.0) %*% t(y.bar - mu.0) / (kappa.0 + n)

    # update Sigma from inverse-Wishart distribution
    sig.iter <- rinvwishart(nu.n, Lambda.n)

    # SWP to calculate conditional parameters
    cond.param <- conditional.parameters(iter.dat) # sweep operator

    ##### I-step

    # for (j in 1:n) {       # row: observations

    if (!is.null(miss.indx)) {

      for (j in 1:length(miss.pos)) {     # col: variables

        miss.p <- miss.pos[j]; miss.in <- miss.indx[, j]
        x <- iter.dat[, miss.p][miss.in == 1]
        x_ <- iter.dat[, -miss.p][miss.in == 1, ]
        mu.x <- mu.iter[miss.p]

        if (length(x) == 1) {
          x <- rnorm(1,
                        mean = mu.x +
                          t(cond.param[miss.p, 2:p]) %*%
                          (x_ - mu.iter[-miss.p]),
                        sd = sqrt(cond.param[miss.p, p + 1]))
        } else {

          for (k in 1:length(x)) {  # impute missing data
            ### decimal places of the imputed values should be the same as the observed values
            x[k] <- rnorm(1,
                          mean = mu.x +
                            t(cond.param[miss.p, 2:p]) %*%
                            (x_[k, ] - mu.iter[-miss.p]),
                          sd = sqrt(cond.param[miss.p, p + 1]))

            # replace the data entry with the imputed data
          }
        }

        iter.dat[, miss.p][miss.in == 1] <- x
      }
    }

    if (!is.null(censor.indx)) {

      if (censor.type == "interval") {

        for (j in 1:length(censor.pos)) {

          censor.p <- censor.pos[j]; censor.in <- censor.indx[, j]
          t. <- iter.dat[, censor.p][censor.in == 1]
          t_ <- iter.dat[, -censor.p][censor.in == 1, ]
          mu.t <- mu.iter[censor.p]

          ll <- censor.val[[1]][, j][censor.in == 1]
          ul <- censor.val[[2]][, j][censor.in == 1]

          if (length(t.) == 1) {

            t. <- rtruncnorm(1,
                             a = ll[k],
                             b = ul[k],
                             mean = mu.t + t(cond.param[censor.p, 2:p]) %*%
                               (t_ - mu.iter[-censor.p]),
                             sd = sqrt(cond.param[censor.p, p + 1]))
          } else {

            for (k in 1:length(t.)){
              t.[k] <- rtruncnorm(1,
                                  a = ll[k],
                                  b = ul[k],
                                  mean = mu.t + t(cond.param[censor.p, 2:p]) %*%
                                    (t_[k, ] - mu.iter[-censor.p]),
                                  sd = sqrt(cond.param[censor.p, p + 1]))
            }
          }
        }
      }
      else if (censor.type == "right") {

        for (j in 1:length(censor.pos)) {

          censor.p <- censor.pos[j]; censor.in <- censor.indx[, j]
          t. <- iter.dat[, censor.p][censor.in == 1]
          t_ <- iter.dat[, -censor.p][censor.in == 1, ]
          mu.t <- mu.iter[censor.p]

          censor.v <- censor.val[, j][censor.in == 1]

          if (length(t.) == 1) {

            t. <- rtruncnorm(1,
                             a = censor.v[k],
                             b = Inf,
                             mean = mu.t + t(cond.param[censor.p, 2:p]) %*%
                               (t_ - mu.iter[-censor.p]),
                             sd = sqrt(cond.param[censor.p, p + 1]))
          } else {

            for (k in 1:length(t.)){
              t.[k] <- rtruncnorm(1,
                                  a = censor.v[k],
                                  b = Inf,
                                  mean = mu.t + t(cond.param[censor.p, 2:p]) %*%
                                    (t_[k, ] - mu.iter[-censor.p]),
                                  sd = sqrt(cond.param[censor.p, p + 1]))
            }
          }
        }
      }
      else if (censor.type == "left") {

        for (j in 1:length(censor.pos)) {

          censor.p <- censor.pos[j]; censor.in <- censor.indx[, j]
          t. <- iter.dat[, censor.p][censor.in == 1]
          t_ <- iter.dat[, -censor.p][censor.in == 1, ]
          mu.t <- mu.iter[censor.p]

          censor.v <- censor.val[, j][censor.in == 1]

          if (length(t.) == 1) {

            t. <- rtruncnorm(1,
                             a = -Inf,
                             b = censor.v[k],
                             mean = mu.t + t(cond.param[censor.p, 2:p]) %*%
                               (t_ - mu.iter[-censor.p]),
                             sd = sqrt(cond.param[censor.p, p + 1]))
          } else {

            for (k in 1:length(t.)){
              t.[k] <- rtruncnorm(1,
                                  a = -Inf,
                                  b = censor.v[k],
                                  mean = mu.t + t(cond.param[censor.p, 2:p]) %*%
                                    (t_[k, ] - mu.iter[-censor.p]),
                                  sd = sqrt(cond.param[censor.p, p + 1]))
            }
          }
        }
      }
      iter.dat[, censor.p][censor.in == 1] <- t.
    }

    # renames
    rownames(sig.iter) <- colnames(sig.iter) <- colnames(dat)

    impute[[i]] <- round(iter.dat, 4)                     # store the imputed dataset for i-th iteration
    Mu.iter[i + 1, ] <- mu.iter                           # store the simulated means from Gibbs Sampler
    Sig.iter[i + 1, ] <- diag(sig.iter)                   # store the simulated variances from Gibbs sampler
    Covmat[[i + 1]] <- sig.iter
    cond.[[i]] <- cond.param

    message(paste(i, "-th iteration!", sep = ""))    # print out the running status

  }
  return(list(
    simulated.mu = Mu.iter,         # simulated mean vector: a vector
    simulated.sig = Sig.iter,       # simulated variance vector: a vector
    simulated.cov = Covmat,         # simulate covariance matrix: a list
    imputed.dat = impute            # simulated data: a list
  ))
}


