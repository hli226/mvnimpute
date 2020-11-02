#' Averaged simulated values plot function
#'
#' Generates average simulated values of the parameter values
#'
#' @param data dataset containing the simulated values.
#' @param iter number of iterations for running multiple imputation.
#' @param x.lab label of the x axis in the generated plot, default is set to "Iteration number".
#' @param y.lab label of the y axis in the generated plot, default is set to "Average of simulated values".
#' @param title title of the generated plot.
#' @param details logical variable to specify whether the average simulated values are returned, default is set to FALSE.
#'
#' @details This function calculates the average simulated values across simulations.
#'
#' @return The plot of averaged values across iterations. If \code{details} = TRUE.
#' A matrix containing the averaged values of all the variables across iterations will be returned.
#'
#' @export
avg.plot <- function(data,
                     iter,
                     x.lab = "Iteration number",
                     y.lab = "Average of simulated values",
                     title = NULL,
                     details = FALSE) {

  # compute the cumulative averages of the parameter values to the current simulation
  avg.param <- matrix(nrow = iter + 1, ncol = ncol(data))
  avg.param[1, ] <- data[1, ]

  for (i in 2:nrow(avg.param)) {

    avg.param[i, ] <- apply(data[1:i, ], 2, mean)

  }

  # transform data from wide to long

  wide.dat <- as.data.frame(cbind(avg.param, 0:iter))
  colnames(wide.dat) <- c(colnames(data), "iter")

  # make sure the id column should be factor variable
  # wide.mean$iter <- as.factor(wide.mean$iter)

  long.dat <- melt(wide.dat, id.vars = "iter",
                   measure.vars = colnames(wide.dat)[colnames(wide.dat) != "iter"],
                   variable.name = "var",
                   value.name = "value")

  print(ggplot(data = long.dat,
         aes(x = as.numeric(iter),
             y = .data$value,
             group = var,
             color = var)) +
    geom_line() +
    xlab(x.lab) +
    ylab(y.lab) +
    ggtitle(title))

  if (details) return(avg.param)

}
