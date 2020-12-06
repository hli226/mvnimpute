####################
## SWEEP Operator ##
####################


swp <- function
(
  dat, # data input
  swp # the index that you want to sweep on
  # NOTE: swp = 0 is to sweep the first row
)
{
  dim.dat <- dim(dat) # dimension of the input dataset
  # calculate the augmented matrix
  p <- dim.dat[2] # number of variables
  n <- dim.dat[1] # number of observations
  aug.mat <- matrix(NA, nrow = p + 1, ncol = p + 1)
  aug.mat[1, 1] <- 1
  for (i in 1:p) {
    aug.mat[1,(i + 1)] <- aug.mat[(i + 1), 1] <- mean(dat[,i])
    # calculate values on the margins, which should be the means of each variable
    for (j in 1:p) {
      aug.mat[(i + 1), (j + 1)] <- aug.mat[(j + 1), (i + 1)] <- sum(dat[, i] * dat[, j]) / n
      # calculate the SS-CP matrix
    }
  }

  aug.mat.dim <- dim(aug.mat)
  rows <- aug.mat.dim[1]
  cols <- aug.mat.dim[2]

  # initiate SWEEP Operator
  swp.mat <- aug.mat
  # create a matrix to store the values after SWEEP Operator
  # this matrix will iterate in the loop

  h.jj <- NA # the diagonal element
  h.ij <- numeric(rows - 1) # the margin vectors
  mat.jk <- matrix(NA, nrow = (rows - 1), ncol = (cols - 1)) # the SS-CP matrix

  for(i in 0:swp) {

    h.jj <- -1 / swp.mat[(i + 1), (i + 1)]
    h.ij <- swp.mat[(i + 1),-(i + 1)] / swp.mat[(i + 1), (i + 1)]
    mat.jk <- swp.mat[-(i + 1), -(i + 1)] - swp.mat[(i + 1),-(i + 1)] %*% t(swp.mat[(i + 1),-(i + 1)]) / swp.mat[(i + 1), (i + 1)]

    swp.mat[i + 1, i + 1] <- h.jj # the swept element
    swp.mat[i + 1, -(i + 1)] <- swp.mat[-(i + 1), i + 1] <- h.ij # the swept row and column
    swp.mat[-(i + 1), - (i+ 1)] <- mat.jk

  }
  return(swp.mat)
}


############################
## Reverse SWEEP Operator ##
############################

Rev.swp <- function(
  dat,
  swp
){

  p <- ncol(dat)

  swp.mat <- swp(dat, p - 1)
  # matrix that has been swept to the greatest extent

  swp.dim <- dim(swp.mat)
  rows <- swp.dim[1]
  cols <- swp.dim[2]

  h.jj <- NA # the diagonal element
  h.ij <- numeric(rows - 1) # the margin vectors
  mat.jk <- matrix(NA, nrow = (rows - 1), ncol = (cols - 1)) # the SS-CP matrix

  for(i in swp:(p - 1)) {

    h.jj <- -1 / swp.mat[(i + 1), (i + 1)]
    h.ij <- -swp.mat[(i + 1),-(i + 1)] / swp.mat[(i + 1), (i + 1)]
    mat.jk <- swp.mat[-(i + 1), -(i + 1)] - swp.mat[(i + 1),-(i + 1)] %*% t(swp.mat[(i + 1),-(i + 1)]) / swp.mat[(i + 1), (i + 1)]

    swp.mat[i + 1, i + 1] <- h.jj # the swept element
    swp.mat[i + 1, -(i + 1)] <- swp.mat[-(i + 1), i + 1] <- h.ij # the swept row and column
    swp.mat[-(i + 1), - (i+ 1)] <- mat.jk

  }
  return(swp.mat)

}

########################################################################################################
## function to calculate the augmented variance-covariance matrix that will be used in SWEEP Operator ##
########################################################################################################

cal.aug.mat <- function (
  dat # data input
)
{
  dim.dat <- dim(dat) # dimension of the input dataset
  # calculate the augmented matrix
  p <- dim.dat[2] # number of variables
  n <- dim.dat[1] # number of observations
  aug.mat <- matrix(NA, nrow = p + 1, ncol = p + 1)
  aug.mat[1, 1] <- 1
  for (i in 1:p) {
    aug.mat[1,(i + 1)] <- aug.mat[(i + 1), 1] <- mean(dat[,i])
    # calculate values on the margins, which should be the means of each variable
    for (j in 1:p) {
      aug.mat[(i + 1), (j + 1)] <- aug.mat[(j + 1), (i + 1)] <- sum(dat[, i] * dat[, j]) / n
      # calculate the SS-CP matrix
    }
  }
  return(aug.mat)
}

#############################
## stepwise SWEEP Operator ##
#############################

step.swp <- function(
  mat, # input matrix to be swept on
  swp # index of row and column to sweep
)
{
  swp.mat <- mat
  i <- swp

  # steps of SWEEP Operator

  h.jj <- -1 / swp.mat[(i + 1), (i + 1)]
  h.ij <- swp.mat[(i + 1),-(i + 1)] / swp.mat[(i + 1), (i + 1)]
  mat.jk <- swp.mat[-(i + 1), -(i + 1)] - swp.mat[(i + 1),-(i + 1)] %*% t(swp.mat[(i + 1),-(i + 1)]) / swp.mat[(i + 1), (i + 1)]

  swp.mat[i + 1, i + 1] <- h.jj # the swept element
  swp.mat[i + 1, -(i + 1)] <- swp.mat[-(i + 1), i + 1] <- h.ij # the swept row and column
  swp.mat[-(i + 1), - (i+ 1)] <- mat.jk

  return(swp.mat)
}

#####################################
## stepwise Reverse SWEEP Operator ##
#####################################


step.rev.swp <- function(
  mat, # inoput matrix to be reversely swept on
  rev.swp # index of row and column to reversely swept
)
{
  rev.swp.mat <- mat
  i <- rev.swp

  # steps of Reverse SWEEP Operator

  h.jj <- -1 / rev.swp.mat[(i + 1), (i + 1)]
  h.ij <- -rev.swp.mat[(i + 1),-(i + 1)] / rev.swp.mat[(i + 1), (i + 1)]
  mat.jk <- rev.swp.mat[-(i + 1), -(i + 1)] - rev.swp.mat[(i + 1),-(i + 1)] %*% t(rev.swp.mat[(i + 1),-(i + 1)]) / rev.swp.mat[(i + 1), (i + 1)]

  rev.swp.mat[i + 1, i + 1] <- h.jj # the swept element
  rev.swp.mat[i + 1, -(i + 1)] <- rev.swp.mat[-(i + 1), i + 1] <- h.ij # the swept row and column
  rev.swp.mat[-(i + 1), - (i+ 1)] <- mat.jk

  return(rev.swp.mat)
}

############################
## conditional parameters ##
############################


conditional.parameters <- function(
  dat
) {

  p <- ncol(dat) # number of variables

  swp.para <- matrix(NA, nrow = p, ncol = p + 1)

  for (i in 1:p) {
    swp.dat <- cbind(dat[, -i], dat[, i]) # i is the index of the variable to be swept on
    swp.para[i, ] <- swp(swp.dat, p - 1)[, p + 1]
  }
  return(swp.para)
}

#############################
## Missing data mechanisms ##
#############################

#### MCAR missing
MCAR.type <- function(
  n,              # number of observations
  p,              # number of variables
  miss.pos,       # variables subject to missing
  miss.percent,   # missing percentages
  censor.pos,     # variables subject to censoring
  censor.percent  # censoring percentages
) {

  #########################################
  ## Generate data type indicator matrix ##
  #########################################

  # create a matrix of the data indicator
  ### 0: missing data
  ### 1: observed data
  ### 2: censored data
  mat <- matrix(1, nrow = n, ncol = p)

  # check if there are variables have both missing and censoring
  common.var <- miss.pos[miss.pos %in% censor.pos]      # the common variables that have both missing and censoring

  ###################################################
  ##  1. a variable has both missing and censoring ##
  ###################################################

  ###### (1) for variable have both missing and censoring

  # only one variable
  if (length(common.var) == 1) {

    j <- which(miss.pos == common.var)
    k <- which(censor.pos == common.var)

    miss.censor.var <- mat[, common.var]
    miss.percent.common <- miss.percent[j]
    censor.percent.common <- censor.percent[k]
    miss.censor.var <- sample(c(0, 1, 2), size = n, replace = TRUE,
                              prob = c(miss.percent.common,
                                       1 - miss.percent.common - censor.percent.common,
                                       censor.percent.common))

    mat[, common.var] <- miss.censor.var

  } else if (length(common.var) > 1) {
    # multiple variables
    for (i in 1:length(common.var)) {

      j <- which(miss.pos == common.var[i])
      k <- which(censor.pos == common.var[i])
      common.var.i <- common.var[i]

      miss.censor.var <- mat[, common.var.i]
      miss.percent.common <- miss.percent[j]
      censor.percent.common <- censor.percent[k]

      miss.censor.var <- sample(c(0, 1, 2), size = n, replace = TRUE,
                                prob = c(miss.percent.common,
                                         1 - miss.percent.common - censor.percent.common,
                                         censor.percent.common))

      mat[, common.var.i] <- miss.censor.var

    }
  }

  ########################################################
  ##  1. for variables have either missing or censoring ##
  ########################################################
  ### (1) missing
  # for variables that have missing other than the common variables

  miss.only <- miss.pos[!miss.pos %in% common.var]

  if (length(miss.only) > 0) {

    ## only one variable
    if (length(miss.only) == 1) {

      j <- which(miss.pos == miss.only)

      miss.only.var <- mat[, miss.only]
      miss.percent.j <- miss.percent[j]

      miss.only.var <- sample(c(0, 1), size = n, replace = TRUE,
                              prob = c(miss.percent.j, 1 - miss.percent.j))

      mat[, miss.only] <- miss.only.var

    } else {

      # multiple variables
      if (length(miss.only) > 1) {

        for (i in 1:length(miss.only)) {

          j <- which(miss.pos == miss.only[i])
          miss.only.var.i <- miss.only[i]

          miss.only.var <- mat[, miss.only.var.i]
          miss.percent.j  <- miss.percent[j]

          miss.only.var <- sample(c(0, 1), size = n, replace = TRUE,
                                  prob = c(miss.percent.j, 1 - miss.percent.j))

          mat[, miss.only.var.i] <- miss.only.var

        }
      }
    }
  }

  ### (2) censoring
  ### for variables that has censoring other than the common variables

  censor.only <- censor.pos[!censor.pos %in% common.var]

  if (length(censor.only > 0)) {

    # only one variable
    if (length(censor.only) == 1) {

      j <- which(censor.pos == censor.only)

      censor.only.var <- mat[, censor.only]
      censor.percent.j <- censor.percent[j]

      censor.only.var <- sample(c(2, 1), size = n, replace = TRUE,
                                prob = c(censor.percent.j, 1 - censor.percent.j))

      mat[, censor.only] <- censor.only.var

    }
    else {

      # multiple variables
      if (length(censor.only) > 1) {

        for (i in 1:length(censor.only)) {

          j <- which(censor.pos == censor.only[i])
          censor.only.var.i <- censor.only[i]

          censor.only.var <- mat[, censor.only.var.i]
          censor.percent.j <- censor.percent[j]

          censor.only.var <- sample(c(2, 1), size = n, replace = TRUE,
                                    prob = c(censor.percent.j, 1 - censor.percent.j))

          mat[, censor.only.var.i] <- censor.only.var

        }
      }
    }
  }

  # rename columns of the resulting matrix
  colnames(mat) <- paste0("y", 1:p)

  return(mat)

}

#### MAR missing
#### data type indicator matrix generation function
MAR.type <- function(
  data,
  miss.pos,       # variables subject to missing
  miss.percent,   # missing percentages
  censor.pos,     # variables subject to censoring
  censor.percent  # censoring percentages
) {

  incomplete.data <- data

  lvalue <- incomplete.data[[1]]; rvalue <- incomplete.data[[2]]

  ## number of observations and variables
  n <- nrow(lvalue); p <- ncol(lvalue)
  data.ind <- matrix(nrow = n, ncol = p)

  # if missing mechnism is MAR, we need to generate observed values first
  var.num <- 1:p
  not.miss <- var.num[!var.num %in% miss.pos]
  obs.pos <- not.miss[!not.miss %in% censor.pos]


  if (length(obs.pos) == 1) {
    # only one observed variable
    incomplete.data[[1]][, obs.pos] <- incomplete.data[[2]][, obs.pos] <- lvalue[, obs.pos]

    # cutpoint of missing data
    # introduce randomness into the missing percentage
    m.percent <- miss.percent * n * runif(1, 0.95, 1.05)
    cutoff <- sort(lvalue[, obs.pos])[m.percent]

    for (i in 1:n) {

      for (j in 1:length(miss.pos)) {
        # for the observed variable has values larger than the cutoff value,
        # the corresponding missing variables will be missing
        if (lvalue[i, obs.pos] < cutoff[j]) {
          incomplete.data[[1]][i, miss.pos[j]] <- -10e10
          incomplete.data[[2]][i, miss.pos[j]] <- 10e10

          # data type indicator
          data.ind[i, miss.pos[j]] <- 0
        }
      }
    }
    data.ind[, obs.pos] <- 1

  }
  else if (length(obs.pos) > 1) {
    # multiple observed variables
    for (i in 1:length(obs.pos)) {
      incomplete.data[[1]][, obs.pos[i]] <- incomplete.data[[2]][, obs.pos[i]] <- lvalue[, obs.pos[i]]
      data.ind[, obs.pos[i]] <- 1
    }

    # missing data depend on the first observed data
    obs.pos.i <- obs.pos[1]

    # cutpoint of missing data
    # introduce randomness into the missing percentage
    m.percent <- miss.percent * n * runif(1, 0.95, 1.05)
    cutoff <- sort(lvalue[, obs.pos.i])[m.percent]

    for (i in 1:n) {

      for (j in 1:length(miss.pos)) {
        # for the observed variable has values larger than the cutoff value,
        # the corresponding missing variables will be missing
        if (lvalue[i, obs.pos.i] < cutoff[j]) {
          incomplete.data[[1]][i, miss.pos[j]] <- -10e10
          incomplete.data[[2]][i, miss.pos[j]] <- 10e10

          data.ind[i, miss.pos[j]] <- 0
        }
      }
    }
  }


  # rename columns of the resulting matrix
  colnames(incomplete.data[[1]]) <- colnames(incomplete.data[[2]]) <- paste0("y", 1:p)

  # return(incomplete.data)
  return(list(data = incomplete.data, ind = data.ind))

}
