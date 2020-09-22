#' Autocorrelation function
#'
#' This function calculates the autocorrelations of all the variables in the dataset
#'
#' @param dat dataset to calculate the autocorrelation
#' @param lag lag of autocorrelation to calulate
#' @param plot  logical variable to specify whether the plot is generated, default is set to FALSE
#' @param title title of the generate autocorrelation plot
#' @param details logical variable to specify whether the autocorrelation values are returned, default is set to TRUE
#'
#' @export
calcu.acf <- function(dat, lag, plot = FALSE, title = NULL, details = TRUE) {

  auto.cor <- matrix(NA, nrow = lag + 1, ncol = ncol(dat))

  for (i in 1:ncol(dat)) {
    # the first row of dat is the initial values
    for (j in 0:lag) auto.cor[j + 1, i] <- round(cor(dat[2:(nrow(dat) - j), i],
                                                     dat[(j + 2):nrow(dat), i]), 3)
  }
  rownames(auto.cor) <- paste("lag ", 0:lag, sep = "")
  colnames(auto.cor) <- colnames(dat)

  if (plot == TRUE) {
    for (i in 1:ncol(auto.cor)) {
      plot(auto.cor[, i], type = "h", ylim = c(-0.5, 1.0),
           xlab = "lag", ylab = "acf",
           main = title[i])
      abline(h = 0)
      abline(h = 0.1, lty = 2, col = "blue")
      abline(h = -0.1, lty = 2, col = "blue")
    }
  }

  if (details == TRUE) return(auto.cor)
}
