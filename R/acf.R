#' Autocorrelation function
#'
#' Calculates the autocorrelations
#'
#' @param dat dataset containing the variables of which autocorrelations are calculated.
#' @param lag lag at which the autocorrelation is calueated, default is set as 50.
#' @param plot  logical variable to specify whether the plot is generated, default is set to FALSE.
#' @param title title of the generated autocorrelation plots.
#' @param details logical variable to specify whether the autocorrelation values are returned, default is set to TRUE.
#'
#' @details This function calculates the autocorrelation values of all the variables on a column by column base. The default value of \code{lag} is set as 50, the maixmum number of lag should not exceed the number of row of the dataset, which reflects the corresponding number of iteration of running the multiple imputation.
#'
#' @return A matrix containing the calculated autocorrelations of all the variables in the dataset. If \code{plot} = TRUE, the autocorrelation plots of all the variables will be drawn.
#'
#' @export
calcu.acf <- function(dat, lag = 50, plot = TRUE, title = NULL, details = FALSE) {

  auto.cor <- matrix(NA, nrow = lag + 1, ncol = ncol(dat))

  for (i in 1:ncol(dat)) {
    # the first row of dat is the initial values
    for (j in 0:lag) auto.cor[j + 1, i] <- round(cor(dat[2:(nrow(dat) - j), i],
                                                     dat[(j + 2):nrow(dat), i]), 3)
  }
  rownames(auto.cor) <- paste("lag ", 0:lag, sep = "")
  colnames(auto.cor) <- colnames(dat)

  if (plot) {
    for (i in 1:ncol(auto.cor)) {
      plot(auto.cor[, i], type = "h", ylim = c(-0.5, 1.0),
           xlab = "lag", ylab = "acf",
           main = title[i])
      abline(h = 0)
    }
  }

  if (details) return(auto.cor)
}
