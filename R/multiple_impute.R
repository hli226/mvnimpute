#' Multiple imputation function
#'
#' Implements the multiple imputation for both missing and censored values
#'
#' @param iter number of rounds for doing multiple imputation
#' @param prior.params list of prior parameter specifications
#' @param initial.values list of starting values
#' @param data dataset to do multiple imputation
#' @param miss.index matrix containing the missing index
#' @param miss.pos vector containing the positions of variables with missing values in the original dataset
#' @param censor.index matrix containing the censoring index
#' @param censor.pos vector containing the position of variables with censored values in the original dataset
#' @param censor.values list containing cutoff values for the censored variables
#' @param censor.type discrete variable specifying the censoring type, it takes "left" for left-censoring, "right" for right-censoring,
#' and "interval" for interval-censoring
#' @param details boolean variable indicating whether the running status is printed in the console. Default set to TRUE
#'
#' @details This function implements the multiple imputation algorithm that concurrently handles missing and censored values. This
#' function requires matrices that contain the respective index of missing or censored values, and vectors that specify the
#' positions of the variables having the missing or censored values in the original dataset. \code{censor.values} includes the cutoff
#' values for the censored values, it should be a matrix which has the same dimension as the censoring index matrix if the data
#' is right- or left-censored, while it should be a list of length 2 if the data is interval-censored. The first element of the
#' list should contain the lower bounds of the censored values, and the upper bounds should be in the second. \code{prior.params}
#' is a list containing the specified parameters of the prior distributions: \code{mu.0} for prior mean; \code{kappa.0} for prior number of measurements,
#' \code{Lambda.0} for the prior scale matrix and \code{nu.0} for the prior degrees of freedom, and \code{initial.values} is a list containing
#' the chosen initial values for the parameters: \code{mu} for the starting values of the mean vector, and \code{sigma} for the starting covariance matrix.
#'
#' @return A list of length 5 including
#'
#' \code{simulated.mu}: The matrix containing the simulated mean values for each variable from each iteration, including the starting values
#'
#' \code{simulated.sig}: The matrix containing the simulated variance values for each variable from each iteration, including the starting values
#'
#' \code{simulated.cov}: A list including the simulated covariance matrix from each iteration
#'
#' \code{imputed.dat}: A list including the imputed data set from each iteration
#'
#' \code{conditional.params}: A list including the conditional parameter values by the sweep operator from each iteration
#'
#' @references
#' Tanner, M., & Wong, W. (1987). The Calculation of Posterior Distributions by Data Augmentation.
#' \emph{Journal of the American Statistical Association}, \bold{82(398)}, 528-540.
#'
#' @export
multiple.impute <- function(
  iter,                           # number of rounds for doing multiple imputation
  prior.params,                   # list of prior parameters
  initial.values,                 # list of initial values
  data,                           # dataset to do imputation (with missing values)
  miss.index = NULL,              # matrix of missing index
  miss.pos = NULL,
  censor.index = NULL,            # matrix of censoring index
  censor.pos = NULL,
  censor.values = NULL,
  censor.type = NULL,
  details = TRUE
)
{

  # input dataset should be a matrix
  if (!is.matrix(data)) stop("The input data matrix should be matrix, please check the type of your data matrix!")

  iter.dat <- data

  # fill.dat. <- initial.impute(dat, miss.pos, censor.pos)

  n <- nrow(iter.dat); p <- ncol(iter.dat)
  # fill the missing values

  # prior parameters
  ### prior specification
  mu.0 <- prior.params$mu.0
  Lambda.0 <- prior.params$Lambda.0
  kappa.0 <- prior.params$kappa.0
  nu.0 <- prior.params$nu.0

  # initial values
  mu.iter <- mu.ini <- initial.values$mu
  sig.iter <- sig.ini <- initial.values$sigma

  # vector and list to store results
  impute <- list()
  Mu.iter <- matrix(nrow = iter + 1, ncol = ncol(iter.dat))
  Sig.iter <- matrix(nrow = iter + 1, ncol = ncol(iter.dat))
  Covmat <- list()
  Mu.iter[1, ] <- mu.ini
  Sig.iter[1, ] <- diag(sig.ini)
  Covmat[[1]] <- sig.ini
  cond. <- list()

  colnames(Mu.iter) <- colnames(Sig.iter) <- colnames(iter.dat)

  # posterior parameters that do not depend on the data
  kappa.n <- kappa.0 + n
  nu.n <- nu.0 + n

  for (i in 1:iter) { ## iteration number

    ###########################################
    ## 1. Use starting values to update data ##
    ###########################################

    # SWP to calculate conditional parameters
    cond.param <- conditional.parameters(iter.dat) # sweep operator

    ##### I-step

    if (!is.null(miss.index)) {

      for (j in 1:length(miss.pos)) {     # col: variables

        miss.p <- miss.pos[j];
        if (!is.null(dim(miss.index))) miss.in <- miss.index[, j] # multivariate missing variables
        else miss.in <- miss.index                                # univariate missing variable

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

    if (!is.null(censor.index)) {

      if (censor.type == "interval") {

        for (j in 1:length(censor.pos)) {

          censor.p <- censor.pos[j]

          if (!is.null(dim(censor.index))) censor.in <- censor.index[, j] # multivariate censored variables
          else censor.in <- censor.index                                  # univariate censored variable

          t. <- iter.dat[, censor.p][censor.in == 1]
          t_ <- iter.dat[, -censor.p][censor.in == 1, ]
          mu.t <- mu.iter[censor.p]

          if (!is.null(dim(censor.index))) {
            ll <- censor.values[[1]][, j][censor.in == 1]                 # multivariate censored variables
            ul <- censor.values[[2]][, j][censor.in == 1]
          } else {
            ll <- censor.values[[1]][censor.in == 1]                      # univariate censored variable
            ul <- censor.values[[2]][censor.in == 1]
          }

          if (length(t.) == 1) {

            t. <- rtruncnorm(1,
                             a = ll,
                             b = ul,
                             mean = mu.t +
                               t(cond.param[censor.p, 2:p]) %*%
                               (t_ - mu.iter[-censor.p]),
                             sd = sqrt(cond.param[censor.p, p + 1]))
          } else {

            for (k in 1:length(t.)){
              t.[k] <- rtruncnorm(1,
                                  a = ll[k],
                                  b = ul[k],
                                  mean = mu.t +
                                    t(cond.param[censor.p, 2:p]) %*%
                                    (t_[k, ] - mu.iter[-censor.p]),
                                  sd = sqrt(cond.param[censor.p, p + 1]))
            }
          }
          iter.dat[, censor.p][censor.in == 1] <- t.
        }
      }

      else if (censor.type == "right") {

        for (j in 1:length(censor.pos)) {

          censor.p <- censor.pos[j]

          if (!is.null(dim(censor.index))) censor.in <- censor.index[, j] # multivariate censored variables
          else censor.in <- censor.index                                  # univariate censored variable

          t. <- iter.dat[, censor.p][censor.in == 1]
          t_ <- iter.dat[, -censor.p][censor.in == 1, ]
          mu.t <- mu.iter[censor.p]

          if (!is.null(dim(censor.index))) {
            censor.v <- censor.values[, j][censor.in == 1]                # multivariate censored variables
          } else {
            censor.v <- censor.values[censor.in == 1]                     # univariate censored variable
          }

          if (length(t.) == 1) {

            t. <- rtruncnorm(1,
                             a = censor.v,
                             b = Inf,
                             mean = mu.t +
                               t(cond.param[censor.p, 2:p]) %*%
                               (t_ - mu.iter[-censor.p]),
                             sd = sqrt(cond.param[censor.p, p + 1]))
          } else {

            for (k in 1:length(t.)){
              t.[k] <- rtruncnorm(1,
                                  a = censor.v[k],
                                  b = Inf,
                                  mean = mu.t +
                                    t(cond.param[censor.p, 2:p]) %*%
                                    (t_[k, ] - mu.iter[-censor.p]),
                                  sd = sqrt(cond.param[censor.p, p + 1]))
            }
          }
          iter.dat[, censor.p][censor.in == 1] <- t.
        }
      }

      else if (censor.type == "left") {

        for (j in 1:length(censor.pos)) {

          censor.p <- censor.pos[j]
          if (!is.null(dim(censor.index))) censor.in <- censor.index[, j] # multivariate censored variables
          else censor.in <- censor.index                                  # univariate censored variable

          t. <- iter.dat[, censor.p][censor.in == 1]
          t_ <- iter.dat[, -censor.p][censor.in == 1, ]
          mu.t <- mu.iter[censor.p]

          if (!is.null(dim(censor.index))) {
            censor.v <- censor.values[, j][censor.in == 1]                # multivariate censored variables
          } else {
            censor.v <- censor.values[censor.in == 1]                     # univariate censored variable
          }

          if (length(t.) == 1) {

            t. <- rtruncnorm(1,
                             a = -Inf,
                             b = censor.v,
                             mean = mu.t +
                               t(cond.param[censor.p, 2:p]) %*%
                               (t_ - mu.iter[-censor.p]),
                             sd = sqrt(cond.param[censor.p, p + 1]))
          } else {

            for (k in 1:length(t.)) {
              t.[k] <- rtruncnorm(1,
                                  a = -Inf,
                                  b = censor.v[k],
                                  mean = mu.t +
                                    t(cond.param[censor.p, 2:p]) %*%
                                    (t_[k, ] - mu.iter[-censor.p]),
                                  sd = sqrt(cond.param[censor.p, p + 1]))
            }

          }
          iter.dat[, censor.p][censor.in == 1] <- t.
        }
      }
    }

    ############################################
    ## 2. Use updated data to update parameters ##
    ############################################

    ### posterior parameters
    y.bar <- apply(iter.dat, 2, mean)
    mu.n <- kappa.0 * mu.0 / (kappa.0 + n) + n  * y.bar / (kappa.0 + n)

    ###### P-step
    # update mu vector from normal distribution condition on Sigma
    mu.iter <- rmvnorm(1, mean = mu.n, sig.iter/kappa.n)

    # posterior parameters for covariance matrix
    S <- apply(iter.dat, 1, "-", mu.iter) %*%
      t(apply(iter.dat, 1, "-", mu.iter))

    Lambda.n <- Lambda.0 + S + kappa.0 * n * (y.bar - mu.0) %*% t(y.bar - mu.0) / (kappa.0 + n)

    # update Sigma from inverse-Wishart distribution
    sig.iter <- rinvwishart(nu.n, Lambda.n)

    # renames
    rownames(sig.iter) <- colnames(sig.iter) <- colnames(iter.dat)

    impute[[i]] <- iter.dat                     # store the imputed dataset for i-th iteration
    Mu.iter[i + 1, ] <- mu.iter                 # store the simulated means from Gibbs sampler
    Sig.iter[i + 1, ] <- diag(sig.iter)         # store the simulated variances from Gibbs sampler
    Covmat[[i + 1]] <- sig.iter                 # store the simulated covariance matrices from Gibbs sampler
    cond.[[i]] <- cond.param                    # store the conditional parameters from SWEEP operator

    if (details) message(paste(i, "-th iteration!", sep = ""))    # print out the running status

  }

  return(list(
    simulated.mu = Mu.iter,         # simulated mean vector: a vector
    simulated.sig = Sig.iter,       # simulated variance vector: a vector
    simulated.cov = Covmat,         # simulate covariance matrix: a list
    imputed.dat = impute,           # simulated data: a list
    conditional.params = cond.      # conditional parameters: a list
  ))

}
