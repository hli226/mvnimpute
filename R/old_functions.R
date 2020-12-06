##### data generation function

data.gen <- function(
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

  # use mvrnorm function from the MASS package to generate multivariate normal data
  dat <- mvrnorm(n, m.vec, sig)

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

    ### censored data
    ### generate censoring times

    censor.ll <- matrix(NA, nrow = n, ncol = censor.num) # lower limit of the censoring interval
    censor.ul <- matrix(NA, nrow = n, ncol = censor.num) # upper limit of the censoring interval
    censor.indx <- matrix(NA, nrow = n, ncol = censor.num) # censoring indices
    # up <- quantile(fully.obs, prob = 0.45)
    # down <- quantile(fully.obs, prob = 0.65)

    for (i in 1:censor.num) {

      censor.p <- censor.pos[i]

      # censored times

      censor.dat <- dat[, censor.p]

      # up <- quantile(fully.obs, prob = 0.45)
      # down <- quantile(fully.obs, prob = 0.65)

      # select n pairs from the censored values as the limits of censoring
      t1 <- runif(n, max(censor.dat) * 0.4, max(censor.dat) * 0.65)
      t2 <- runif(n, max(censor.dat) * 0.7, max(censor.dat) * 0.88)

      # t1 <- rexp(n, 1 / up)
      # t2 <- rexp(n , 1 / down)
      censor.ll[, i] <- t1
      censor.ul[, i] <- t2

      # generating censoring indices

      # censor.dat <- dat[, censor.p]
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
      miss.prob <- sample(seq(0.1, 0.5, by = 0.001), 1, replace = TRUE)
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
    # up <- quantile(fully.obs, prob = 0.45)
    # down <- quantile(fully.obs, prob = 0.65)

    for (i in 1:censor.num) {
      censor.p <- censor.pos[i]

      censor.dat <- dat[, censor.p]

      # t1 <- rexp(n, 1 / up)
      # t2 <- rexp(n , 1 / down)

      t1 <- runif(n, max(censor.dat) * 0.4, max(censor.dat) * 0.65)
      t2 <- runif(n, max(censor.dat) * 0.7, max(censor.dat) * 0.88)

      censor.ll[, i] <- t1
      censor.ul[, i] <- t2

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

##### single imputation function

initial.impute <- function(data,
                           miss.index = NULL, miss.pos = NULL,
                           censor.index = NULL, censor.pos = NULL
)
{
  miss.dat <- data

  # fill-in missing data

  if (!is.null(miss.index)){
    for (i in 1:length(miss.pos)) {

      miss <- miss.pos[i]
      if (!is.null(dim(miss.index))) miss.in <- miss.index[, i]
      else miss.in <- miss.index

      x <- miss.dat[, miss][miss.in == 1]
      miss.mean <- calcu.param(miss.dat)$CC.mean[miss]
      miss.var <- calcu.param(miss.dat)$CC.var[miss]

      for (j in 1:length(x)) {
        x[j] <- rnorm(1, miss.mean, sqrt(miss.var))
        # generate normal random variable using the complete cases mean and variance
      }

      message(paste("Impute missing data column ", i, sep = ""))

      miss.dat[, miss][miss.in == 1] <- x

    }
  }


  # fill-in censored data
  if (!is.null(censor.index)) {
    for (i in 1:length(censor.pos)) {

      censor <- censor.pos[i]
      if (!is.null(dim(censor.index))) censor.in <- censor.index[, i]
      else censor.in <- censor.index

      t. <- miss.dat[, censor][censor.in == 1]
      censor.mean <- calcu.param(miss.dat)$CC.mean[censor]
      censor.var <- calcu.param(miss.dat)$CC.var[censor]

      for (j in 1:length(t.)) {

        t.[j] <- rnorm(1,
                       mean = censor.mean, sqrt(censor.var))
        # generate random variable uniformly in the interval
      }

      message(paste("Impute censored data column ", i, sep = ""))

      miss.dat[, censor][censor.in == 1] <- t.

    }
  }

  return(as.matrix(miss.dat))
}

##### multiple imputation function

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
    mu.iter <- mvrnorm(1, mu.n, sig.iter/kappa.n)

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


##### CC and AC parameter calculations function
calcu.param <- function(data) {

  if (is.null(dim(data))) stop("Need a multidimensional matrix!")

  ## CC parameters
  complete.dat <- data[complete.cases(data), ]
  CC.mean <- apply(complete.dat, 2, mean)
  CC.var <- apply(complete.dat, 2, var)
  CC.cov <- cov(complete.dat)

  ## AC parameters
  AC.mean <- numeric(ncol(data)); AC.var <- numeric(ncol(data))
  for (i in 1:ncol(data)) {
    available.dat <- na.omit(data[, i])
    AC.mean[i] <- mean(available.dat)
    AC.var[i] <- var(available.dat)
  }

  return(list(
    AC.mean = AC.mean,
    AC.var = AC.var,
    CC.cov = CC.cov,
    CC.mean = CC.mean,
    CC.var = CC.var))
}

##### visualization plot

visual.plot <- function(
  data,
  miss.index = NULL,
  censor.index = NULL,
  title = NULL) {

  `%notin%` <- Negate(`%in%`)

  # 1. data only has missing values
  if (!is.null(miss.index) & is.null(censor.index)) {
    # either missing or censoring index is null matrix
    plot.dat <- data
    miss.indx <- as.data.frame(miss.index)

    missing.vals <- miss.indx %>%
      gather(key = "key", value = "val") %>%
      mutate(isna = (.data$val == 1)) %>%
      group_by(.data$key) %>%
      mutate(total = n()) %>%
      group_by(.data$key, .data$total, .data$isna) %>%
      summarise(num.isna = n()) %>%
      mutate(pct = .data$num.isna / .data$total * 100)

    no.missing <- colnames(plot.dat)[colnames(plot.dat) %notin% unique(missing.vals$key)]

    if (length(no.missing) > 0) {

      no.missingfalse <- data.frame(
        key = no.missing,
        total = nrow(plot.dat),
        isna = FALSE,
        num.isna = nrow(plot.dat),
        pct = 100
      )

      no.missingtrue <- data.frame(
        key = no.missing,
        total = nrow(plot.dat),
        isna = TRUE,
        num.isna = 0,
        pct = 0
      )

      missing.vals <- rbind(missing.vals, no.missingfalse, no.missingtrue)

    }

    levels <- (missing.vals %>% filter(.data$isna == TRUE) %>% arrange(desc(.data$pct)))$key

    missing.vals %>%
      ggplot() +
      geom_bar(aes(x = reorder(.data$key, desc(.data$pct)),
                   y = .data$pct, fill = .data$isna),
               stat = "identity") +
      scale_x_discrete(limits = levels,
                       labels = paste(
                         levels,
                         " (",
                         round((missing.vals %>% filter(.data$isna == TRUE) %>% arrange(desc(.data$pct)))$pct, 2),
                         "%)", sep = ""
                       )) +
      scale_fill_manual(name = "",
                        values = c("steelblue", "tomato3"),
                        labels = c("Observed", "Missing")) +
      coord_flip() +
      labs(
        x = "Variable",
        y = "Percentage") +
      ggtitle(title)

  }

  # 2. data only has censored values
  else if (is.null(miss.index) & !is.null(censor.index)) {

    plot.dat <- data
    censor.indx <- as.data.frame(censor.index)

    censoring.vals <- censor.indx %>%
      gather(key = "key", value = "val") %>%
      mutate(iscensor = (.data$val == 1)) %>%
      group_by(.data$key) %>%
      mutate(total = n()) %>%
      group_by(.data$key, .data$total, .data$iscensor) %>%
      summarise(num.iscensor = n()) %>%
      mutate(pct = .data$num.iscensor / .data$total * 100)

    no.censoring <- colnames(plot.dat)[colnames(plot.dat) %notin% unique(censoring.vals$key)]

    if (length(no.censoring) > 0) {

      no.censoringfalse <- data.frame(
        key = no.censoring,
        total = nrow(plot.dat),
        iscensor = FALSE,
        num.iscensor = nrow(plot.dat),
        pct = 100
      )

      no.censoringtrue <- data.frame(
        key = no.censoring,
        total = nrow(plot.dat),
        iscensor = TRUE,
        num.iscensor = 0,
        pct = 0
      )

      censoring.vals <- rbind(censoring.vals, no.censoringfalse, no.censoringtrue)

    }

    levels <- (censoring.vals %>% filter(.data$iscensor == TRUE) %>% arrange(desc(.data$pct)))$key

    censoring.vals %>%
      ggplot() +
      geom_bar(aes(x = reorder(.data$key, desc(.data$pct)),
                   y = .data$pct, fill = .data$iscensor),
               stat = "identity") +
      scale_x_discrete(limits = levels,
                       labels = paste(
                         levels,
                         " (",
                         round((censoring.vals %>% filter(.data$iscensor == TRUE) %>% arrange(desc(.data$pct)))$pct, 2),
                         "%)", sep = ""
                       )) +
      scale_fill_manual(name = "",
                        values = c("steelblue", "lightgreen"),
                        labels = c("Observed", "Censored")) +
      coord_flip() +
      labs(
        x = "Variable",
        y = "Percentage") +
      ggtitle(title)

  }

  else if (!is.null(miss.index) & !is.null(censor.index)) {

    plot.dat <- as.data.frame(data)
    miss.indx <- as.data.frame(miss.index)
    censor.indx <- as.data.frame(censor.index)

    miss.name <- colnames(miss.indx)
    censor.name <- colnames(censor.indx)

    common.name  <- NULL

    for (i in 1:length(miss.name)) {

      if (miss.name[i] %in% censor.name) common.name <- c(common.name, miss.name[i])

    }

    # (1) if all the variable have either missing or censored values
    if (is.null(common.name)) {

      missing.vals <- miss.indx %>%
        gather(key = "key", value = "val") %>%
        mutate(isna = (.data$val == 1)) %>%
        group_by(.data$key) %>%
        mutate(total = n()) %>%
        group_by(.data$key, .data$total, .data$isna) %>%
        summarise(num.isna = n()) %>%
        mutate(pct = .data$num.isna / .data$total * 100)

      censoring.vals <- censor.indx %>%
        gather(key = "key", value = "val") %>%
        mutate(iscensor = (.data$val == 1)) %>%
        group_by(.data$key) %>%
        mutate(total = n()) %>%
        group_by(.data$key, .data$total, .data$iscensor) %>%
        summarise(num.iscensor = n()) %>%
        mutate(pct = .data$num.iscensor / .data$total * 100)

      missing.vals <- data.frame(missing.vals, type = ifelse(missing.vals$isna == TRUE, "0", "1"))
      missing.vals <- missing.vals[, -3]
      colnames(missing.vals) <- c("key", "total", "num", "pct", "type")

      censoring.vals <- censoring.vals[order(censoring.vals[, 1]), ]

      censoring.vals <- data.frame(censoring.vals, type = ifelse(censoring.vals$iscensor == TRUE, "2", "1"))
      censoring.vals <- censoring.vals[, -3]
      colnames(censoring.vals) <- c("key", "total", "num", "pct", "type")

      observing <- colnames(plot.dat)[colnames(plot.dat) %notin% c(unique(missing.vals$key), unique(censoring.vals$key))]

      if (length(observing) > 0) {

        observing.vals <- data.frame(
          key = observing,
          total = nrow(plot.dat),
          num = rep(nrow(plot.dat), length(observing)),
          pct = 100,
          type = 1
        )
      }
      complete.dat <- rbind(missing.vals, censoring.vals, observing.vals)

    }

    # if some of the variables have both missing and censored values
    else if (!is.null(common.name)) {

      missing.vals <- miss.indx %>%
        gather(key = "key", value = "val") %>%
        mutate(isna = (.data$val == 1)) %>%
        group_by(.data$key) %>%
        mutate(total = n()) %>%
        group_by(.data$key, .data$total, .data$isna) %>%
        summarise(num.isna = n()) %>%
        mutate(pct = .data$num.isna / .data$total * 100) %>%
        filter(.data$isna == TRUE)

      missing.other <- miss.indx %>%
        gather(key = "key", value = "val") %>%
        mutate(isna = (.data$val == 1)) %>%
        group_by(.data$key) %>%
        mutate(total = n()) %>%
        group_by(.data$key, .data$total, .data$isna) %>%
        summarise(num.isna = n()) %>%
        mutate(pct = .data$num.isna / .data$total * 100) %>%
        filter(.data$isna == FALSE)

      no.missing <- missing.other$key[missing.other$key %notin% missing.vals$key]

      if (length(no.missing) > 0) {

        missing.0 <- data.frame(
          key = no.missing,
          total = nrow(plot.dat),
          isna = rep(TRUE, length(no.missing)),
          num.isna = rep(0, length(no.missing)),
          pct = rep(0, length(no.missing))
        )

        missing.vals <- rbind(missing.vals, missing.0)
      }

      censoring.vals <- censor.indx %>%
        gather(key = "key", value = "val") %>%
        mutate(iscensor = (.data$val == 1)) %>%
        group_by(.data$key) %>%
        mutate(total = n()) %>%
        group_by(.data$key, .data$total, .data$iscensor) %>%
        summarise(num.iscensor = n()) %>%
        mutate(pct = .data$num.iscensor / .data$total * 100) %>%
        filter(.data$iscensor == TRUE)

      censoring.other <- censor.indx %>%
        gather(key = "key", value = "val") %>%
        mutate(iscensor = (.data$val == 1)) %>%
        group_by(.data$key) %>%
        mutate(total = n()) %>%
        group_by(.data$key, .data$total, .data$iscensor) %>%
        summarise(num.iscensor = n()) %>%
        mutate(pct = .data$num.iscensor / .data$total * 100) %>%
        filter(.data$iscensor == FALSE)

      no.censoring <- censoring.other$key[censoring.other$key %notin% censoring.vals$key]

      if (length(no.censoring) > 0) {

        censoring.0 <- data.frame(
          key = no.censoring,
          total = nrow(plot.dat),
          iscensor = rep(TRUE, length(no.censoring)),
          num.iscensor = rep(0, length(no.censoring)),
          pct = rep(0, length(no.censoring))
        )

        censoring.vals <- rbind(censoring.vals, censoring.0)

      }

      missing.vals <- data.frame(missing.vals, type = rep("0", nrow(missing.vals)))
      missing.vals <- missing.vals[, -3]
      colnames(missing.vals) <- c("key", "total", "num", "pct", "type")
      missing.vals <- missing.vals[order(missing.vals$key), ]

      censoring.vals <- data.frame(censoring.vals, type = ifelse(censoring.vals$iscensor == TRUE, "2", "1"))
      censoring.vals <- censoring.vals[, -3]
      colnames(censoring.vals) <- c("key", "total", "num", "pct", "type")
      censoring.vals <- censoring.vals[order(censoring.vals$key), ]

      observing.vals <- data.frame(
        key = missing.vals$key,
        total = rep(nrow(plot.dat), length(missing.vals$key)),
        num = nrow(plot.dat) - missing.vals$num - censoring.vals$num,
        pct = (nrow(plot.dat) - missing.vals$num - censoring.vals$num) / nrow(plot.dat),
        type = rep(1, length(missing.vals$key))
      )

      complete.dat <- rbind(missing.vals,
                            censoring.vals,
                            observing.vals)

      complete.dat$pct <- complete.dat$num / complete.dat$total * 100

    }

    complete.dat <- complete.dat %>%
      arrange(.data$key, .data$type)

    complete.dat %>%
      ggplot() +
      geom_bar(aes(x = .data$key,
                   y = as.numeric(.data$pct),
                   fill = .data$type), stat = "identity") +
      scale_fill_manual(name = "",
                        values = c("tomato3", "steelblue", "lightgreen"),
                        labels = c("Missing", "Observed", "Censored")) +
      coord_flip() +
      xlab("Variable") +
      ylab("Percentage") +
      ggtitle(title)

  }

}

