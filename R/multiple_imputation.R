#' Multiple imputation function
#'
#' Function to do multiple imputations
#'
#' @param iter Number of iterations
#' @param prior.param A list containing the parameter values of the prior distribution
#' @param ini.vals A list containing the initial values
#' @param dat Dataset to do multiple imputation
#' @param missing.indx A matrix containing the missing indices
#' @param censoring.indx A matrix containing the censoring indices
#' @param censoring.val A matrix containing the censoring information
#' @param censoring.type Categorical variable specifying whether the data is right censored, left censored, or interval censored
#'
#' @export
sim.func <- function(
  iter,                       # number of iterations to run
  prior.param,                # list of prior parameters
  ini.vals,                   # list of initial values
  dat,                        # dataset to do imputation (with missing values)
  censoring.val,              # matrix of censoring values
  missing.indx = NULL,        # matrix of missing index
  censoring.indx = NULL,      # matrix of censoring index
  censoring.type = "left"     # specify the type of censoring, default is set to interval censoring
)
{
  n <- nrow(dat); p <- ncol(dat)
  # fill the missing values
  iter.dat <- fill.dat(dat)

  # prior parameters
  ### prior specification
  mu.0 <- prior.param[[1]]
  V <- prior.param[[2]]

  kappa.0 <- prior.param[[3]]
  nu.0 <- prior.param[[4]]; Lambda.0 <- prior.param[[5]]
  rownames(V) <- colnames(V) <- NULL

  # initial values
  mu.iter <- mu.ini <- ini.vals[[1]]
  sig.iter <- sig.ini <- ini.vals[[2]]

  # vector and list to store results
  impute <- list()
  Mu.iter <- matrix(nrow = iter + 1, ncol = ncol(dat))
  Sig.iter <- matrix(nrow = iter + 1, ncol = ncol(dat))
  Covmat <- list()
  Mu.iter[1, ] <- mu.ini
  Sig.iter[1, ] <- diag(sig.ini)
  Covmat[[1]] <- sig.ini

  colnames(Mu.iter) <- colnames(Sig.iter) <- colnames(dat)

  # posterior parameters that do not depend on the data
  kappa.n <- kappa.0 + n
  nu.n <- nu.0 + n

  # simulation starting time
  start <- Sys.time()

  # # setting seed
  # set.seed(seed)

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

      for (k in 1:p) {     # column: variables

        if (!is.null(missing.indx)) {

          if (missing.indx[j, k] == 1) {  # impute missing data
            ### decimal places of the imputed values should be the same as the observed values
            miss.dat <- rnorm(1,
                              mean = mu.iter[k] +
                                t(cond.param[k, 2:p]) %*% (iter.dat[j, -k] - mu.iter[-k]),
                              sd = sqrt(cond.param[k, p+1]))

            # replace the data entry with the imputed data
            iter.dat[j, k] <- miss.dat
          }
        }

        if (!is.null(censoring.indx)) {

          if (censoring.type == "left") {

            if (censoring.indx[j, k] == 1) { # impute censored data

              censor.dat <- rtruncnorm(1,
                                       a = -Inf,
                                       b = censoring.val[j, k],
                                       mean = mu.iter[k] +
                                         t(cond.param[k, 2:p]) %*% (iter.dat[j, -k] - mu.iter[-k]),
                                       sd = sqrt(cond.param[k, p+1]))
              # replace the data entry with the imputed data

              iter.dat[j, k] <- censor.dat
              # what are the indices of censoring data
            }
          }
        }
      }
    }
    # renames
    rownames(sig.iter) <- colnames(sig.iter) <- colnames(iter.dat) <- colnames(dat)

    impute[[i]] <- round(iter.dat, 4)                 # store the imputed dataset for i-th iteration
    Mu.iter[i + 1, ] <- mu.iter                           # store the simulated means from Gibbs Sampler
    Sig.iter[i + 1, ] <- diag(sig.iter)                   # store the simulated variances from Gibbs sampler
    Covmat[[i + 1]] <- sig.iter

    message(paste(i, "-th iteration!", sep = ""))    # print out the running status
  }
  # simulation ending time
  end <- Sys.time()

  # simulation running time
  dur <- paste("The program takes ", end - start, " to run.", sep = "")

  return(list(
    simulated.mu = Mu.iter,         # simulated mean vector: a vector
    simulated.sig = Sig.iter, # simulated variance vector: a vector
    simulated.cov = Covmat,        # simulate covariance matrix: a list
    imputed.dat = impute,           # simulated data: a list
    running.time = dur              # running time
  ))
}
