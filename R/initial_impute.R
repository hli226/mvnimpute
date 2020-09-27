#' Initial missing data fill-in fucntion
#'
#' Fills in data with missing and censored values for multiple imputation
#'
#' @param dat dataset with missing and censored values
#' @param method method to be used for initial filled-in values
#' @param miss.pos Position of variables with missing values in the original dataset
#' @param censor.pos Position of variables with censored values in the original dataset
#'
#' @details The missing values will initially be filled in by some single imputation methods. Currently, it only supports generating
#' random values from the normal distribution with the mean and variance as the complete-case mean and variance, respectively.
#'
#' @return a complete dataset with missing values filled in
#' @export
initial.impute <- function(dat, miss.pos, censor.pos, method = "normal")
{

  fill.dat <- dat

  if (method == "normal") {

    m <- calcu.param(dat)$CC.mean; v <- calcu.param(dat)$CC.var

    for (i in 1:nrow(dat)) {

      for (j in 1:length(censor.pos)) {
        fill.dat[i, censor.pos[j]] <- rnorm(1, m[censor.pos[j]], sqrt(v[censor.pos[j]]))
      }

      for (k in 1:length(miss.pos)) {
        fill.dat[i, miss.pos[k]] <- rnorm(1, m[miss.pos[k]], sqrt(v[miss.pos[k]]))
      }
    }

  }
  return(fill.dat)
}
