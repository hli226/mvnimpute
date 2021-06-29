#' Multiple imputation function
#'
#' Implements the multiple imputation for both missing and censored values
#'
#' @param data a list of data containing the information for the missing and censored values
#' @param prior.params list of prior parameter specifications
#' @param starting.values list of starting values
#' @param iter number of rounds for doing multiple imputation
#' @param details boolean variable indicating whether the running status is printed in the console. Default is set to TRUE
#'
#' @references
#' Tanner, M., & Wong, W. (1987). The Calculation of Posterior Distributions by Data Augmentation.
#' \emph{Journal of the American Statistical Association}, \bold{82(398)}, 528-540.
#'
#' @export
multiple.imputation <- function(
  data,          # the list that contains the censored values
  prior.params,  # prior specifications
  starting.values,# starting values
  iter,          # iterations of Gibbs sampler
  details = TRUE # boolean variable to print out running status
) {

  ### single imputation to make up incomplete data
  iter.data <- single.imputation(data)

  lvalue <- data[[1]]; rvalue <- data[[2]]
  # number of observations and variables
  n <- nrow(lvalue); p <- ncol(lvalue)

  ###########################################
  ### prior specification
  mu.0 <- prior.params$mu.0
  Lambda.0 <- prior.params$Lambda.0
  kappa.0 <- prior.params$kappa.0
  nu.0 <- prior.params$nu.0

  # starting values
  mu.iter <- mu.ini <- starting.values$mu
  sig.iter <- sig.ini <- starting.values$sigma

  # vector and list to store results
  impute <- list()
  Mu.iter <- matrix(nrow = iter + 1, ncol = ncol(iter.data))
  Sig.iter <- matrix(nrow = iter + 1, ncol = ncol(iter.data))
  Covmat <- list()

  ### set the first set of values as the staring values
  Mu.iter[1, ] <- mu.ini
  Sig.iter[1, ] <- diag(sig.ini)
  Covmat[[1]] <- sig.ini
  cond. <- list()

  # rename
  colnames(Mu.iter) <- colnames(Sig.iter) <- colnames(lvalue)

  # posterior parameters that do not depend on the data
  kappa.n <- kappa.0 + n
  nu.n <- nu.0 + n

  ## Gibbs iteration
  for (i in 1:iter) {

    ###########################################
    ## 1. Use starting values to update data ##
    ###########################################

    # SWP to calculate conditional parameters
    cond.param <- conditional.parameters(iter.data) # sweep operator

    ##### I-step

    for (j in 1:n) {

      for (k in 1:p) {
        # (1) missing values
        if (lvalue[j, k] == -10e10 & rvalue[j, k] == 10e10) {
          iter.data[j, k] <- rnorm(1,
                                   mean = mu.iter[k] + t(cond.param[k, 2:p]) %*%
                                     (iter.data[j, -k] - mu.iter[-k]),
                                   sd = sqrt(cond.param[k, p + 1]))
        }
        else if (lvalue[j, k] != rvalue[j, k] &
                 (lvalue[j, k] > -10e10 & rvalue[j, k] < 10e10)) {
          # (2) censoring values
          # (2.1) interval censoring
          iter.data[j, k] <- rtruncnorm(1,
                                        a = lvalue[j, k],
                                        b = rvalue[j, k],
                                        mean = mu.iter[k] + t(cond.param[k, 2:p]) %*%
                                          (iter.data[j, -k] - mu.iter[-k]),
                                        sd = sqrt(cond.param[k, p + 1]))
        }
        else if (lvalue[j, k] != rvalue[j, k] &
                 (lvalue[j, k] > -10e10 & rvalue[j, k] == 10e10))  {
          # (2.2) right censoring
          iter.data[j, k] <- rtruncnorm(1,
                                        a = lvalue[j, k],
                                        mean = mu.iter[k] + t(cond.param[k, 2:p]) %*%
                                          (iter.data[j, -k] - mu.iter[-k]),
                                        sd = sqrt(cond.param[k, p + 1]))
        }
        else if (lvalue[j, k] != rvalue[j, k] &
                 (lvalue[j, k] == -10e10 & rvalue[j, k] < 10e10)) {
          # (2.3) left censoring
          iter.data[j, k] <- rtruncnorm(1,
                                        b = rvalue[j, k],
                                        mean = mu.iter[k] + t(cond.param[k, 2:p]) %*%
                                          (iter.data[j, -k] - mu.iter[-k]),
                                        sd = sqrt(cond.param[k, p + 1]))
        }
      }
    }

    ##############################################
    ## 2. Use updated data to update parameters ##
    ##############################################

    ### posterior parameters
    y.bar <- apply(iter.data, 2, mean)
    mu.n <- kappa.0 * mu.0 / (kappa.0 + n) + n  * y.bar / (kappa.0 + n)

    ###### P-step
    # update mu vector from normal distribution condition on Sigma
    mu.iter <- mvrnorm(1, mu.n, sig.iter/kappa.n)

    # posterior parameters for covariance matrix
    S <- apply(iter.data, 1, "-", mu.iter) %*%
      t(apply(iter.data, 1, "-", mu.iter))

    Lambda.n <- Lambda.0 + S + kappa.0 * n * (y.bar - mu.0) %*% t(y.bar - mu.0) / (kappa.0 + n)

    # update Sigma from inverse-Wishart distribution
    sig.iter <- rinvwishart(nu.n, Lambda.n)

    # rename
    rownames(sig.iter) <- colnames(sig.iter) <- colnames(lvalue)

    impute[[i]] <- iter.data                                      # store the imputed dataset for i-th iteration
    Mu.iter[i + 1, ] <- mu.iter                                   # store the simulated means from Gibbs sampler
    Sig.iter[i + 1, ] <- diag(sig.iter)                           # store the simulated variances from Gibbs sampler
    Covmat[[i + 1]] <- sig.iter                                   # store the simulated covariance matrices from Gibbs sampler
    cond.[[i]] <- cond.param                                      # store the conditional parameters from SWEEP operator
    if (details) message(paste(i, "-th iteration!", sep = ""))    # print out the running status
  }

  return(list(
    simulated.mu = Mu.iter,         # simulated mean vector: a vector
    simulated.sig = Sig.iter,       # simulated variance vector: a vector
    simulated.cov = Covmat,         # simulate covariance matrix: a list
    impute.data = impute,           # simulated data: a list
    conditional.params = cond.      # conditional parameters: a list
  ))

}
