#' Multiple imputation function
#'
#' This function implements the multiple imputation of the missing and censored values
#'
#' @param iter Iteration number
#' @param prior.param List of prior parameter values
#' @param ini.vals List of initial values
#' @param dat Dataset to do multiple imputation
#' @param miss.indx Matrix of missing data index
#' @param miss.pos Position of variables with missing values in the original dataset
#' @param censor.indx Matrix of censoring data index
#' @param censor.pos Position of variables with censored values in the orignal dataset
#' @param censor.val List containing cutoff values for the censored variables
#' @param censor.type Logical variable specifying the censoring type
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

  fill.dat. <- initial.impute(dat, miss.pos, censor.pos)

  n <- nrow(fill.dat.); p <- ncol(fill.dat.)
  # fill the missing values

  # prior parameters
  ### prior specification
  mu.0 <- prior.param[[1]]
  V.0 <- prior.param[[2]]

  kappa.0 <- prior.param[[3]]
  nu.0 <- prior.param[[4]]; Lambda.0 <- prior.param[[5]]

  # initial values
  mu.iter <- mu.ini <- ini.vals[[1]]
  sig.iter <- sig.ini <- ini.vals[[2]]

  # vector and list to store results
  impute <- list()
  Mu.iter <- matrix(nrow = iter + 1, ncol = ncol(fill.dat.))
  Sig.iter <- matrix(nrow = iter + 1, ncol = ncol(fill.dat.))
  Covmat <- list()
  Mu.iter[1, ] <- mu.ini
  Sig.iter[1, ] <- diag(sig.ini)
  Covmat[[1]] <- sig.ini
  cond. <- list()

  colnames(Mu.iter) <- colnames(Sig.iter) <- colnames(dat)

  # posterior parameters that do not depend on the data
  kappa.n <- kappa.0 + n
  nu.n <- nu.0 + n

  iter.dat <- fill.dat.

  for (i in 1:iter) { ## iteration number

    ### posterior parameters
    y.bar <- apply(iter.dat, 2, mean)
    mu.n <- kappa.0 * mu.0 / (kappa.0 + n) + n  * y.bar / (kappa.0 + n)

    # P-step
    # update mu vector from normal distribution condition on Sigma
    mu.iter <- rmvnorm(1, mean = mu.n, sig.iter/kappa.n)

    # posterior parameters for variance-covariance matrix
    S <- apply(iter.dat, 1, "-", mu.iter) %*%
      t(apply(iter.dat, 1, "-", mu.iter))

    Lambda.n <- Lambda.0 + S + kappa.0 * n * (y.bar - mu.0) %*% t(y.bar - mu.0) / (kappa.0 + n)

    # update Sigma from inverse-Wishart distribution
    sig.iter <- rinvwishart(nu.n, Lambda.n)

    # SWP to calculate conditional parameters
    cond.param <- conditional.parameters(iter.dat) # sweep operator

    ##### I-step

    for (j in 1:n) {       # row: observations

      if (!is.null(miss.indx)) {

        for (k in 1:ncol(miss.indx)) {     # col: variables

          if (miss.indx[j, k] == 1) {  # impute missing data
            ### decimal places of the imputed values should be the same as the observed values
            miss.dat <- rnorm(1,
                              mean = mu.iter[miss.pos][k] +
                                t(cond.param[miss.pos[k], 2:p]) %*%
                                (iter.dat[j, -miss.pos[k]] - mu.iter[-miss.pos[k]]),
                              sd = sqrt(cond.param[miss.pos[k], p + 1]))

            # replace the data entry with the imputed data
            iter.dat[j, miss.pos[k]] <- miss.dat
          }
        }
      }

      if (!is.null(censor.indx)) {

        ## interval censoring

        for (l in 1:ncol(censor.indx)) {

          if (censor.type == "interval") {
            if (censor.indx[j, l] == 1) { # impute censored data

              censor.dat <- rtruncnorm(1,
                                       a = censor.val[[1]][j, l],
                                       b = censor.val[[2]][j, l],
                                       mean = mu.iter[censor.pos][l] +
                                         t(cond.param[censor.pos[l], 2:p]) %*%
                                         (iter.dat[j, -censor.pos[l]] - mu.iter[-censor.pos[l]]),
                                       sd = sqrt(cond.param[censor.pos[l], p + 1]))

              # replace the data entry with the imputed data
              iter.dat[j, censor.pos[l]] <- censor.dat
            }
          }
          # left censoring

          else if (censor.type == "left") {

          if (censor.indx[j, l] == 1) { # impute censored data

            censor.dat <- rtruncnorm(1,
                                     a = -Inf,
                                     b = censor.val[j, l],
                                     mean = mu.iter[censor.pos][l] +
                                       t(cond.param[censor.pos[l], 2:p]) %*%
                                       (iter.dat[j, -censor.pos[l]] - mu.iter[-censor.pos[l]]),
                                     sd = sqrt(cond.param[censor.pos[l], p + 1]))

            # replace the data entry with the imputed data
            iter.dat[j, censor.pos[l]] <- censor.dat

          }

          }
          ## right censoring
          if (censor.type == "right") {
          if (censor.indx[j, l] == 1) { # impute censored data
            censor.dat <- rtruncnorm(1,
                                     a = censor.val[j, l],
                                     b = Inf,
                                     mean = mu.iter[censor.pos][l] +
                                       t(cond.param[censor.pos[l], 2:p]) %*%
                                       (iter.dat[j, -censor.pos[l]] - mu.iter[-censor.pos[l]]),
                                     sd = sqrt(cond.param[censor.pos[l], p + 1]))

            # replace the data entry with the imputed data
            iter.dat[j, censor.pos[l]] <- censor.dat
           }
          }
        }
      }
    }
    # renames
    rownames(sig.iter) <- colnames(sig.iter) <- colnames(dat)

    impute[[i]] <- round(iter.dat, 4)                 # store the imputed dataset for i-th iteration
    Mu.iter[i + 1, ] <- mu.iter                           # store the simulated means from Gibbs Sampler
    Sig.iter[i + 1, ] <- diag(sig.iter)                   # store the simulated variances from Gibbs sampler
    Covmat[[i + 1]] <- sig.iter
    cond.[[i]] <- cond.param

    message(paste(i, "-th iteration!", sep = ""))    # print out the running status

  }


  return(list(
    simulated.mu = Mu.iter,         # simulated mean vector: a vector
    simulated.sig = Sig.iter, # simulated variance vector: a vector
    simulated.cov = Covmat,        # simulate covariance matrix: a list
    imputed.dat = impute        # simulated data: a list
  ))

  # return(cond.)

}

