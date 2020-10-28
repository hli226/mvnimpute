#' Convergence plot function
#'
#' Generates the convergence plots for the parameter values
#'
#' @param data dataset containing the simulated values.
#' @param iter number of iterations for running multiple imputation.
#' @param x.lab label of the x axis in the generated plot, default is set to "Iteration number".
#' @param y.lab label of the y axis in the generated plot, default is set to "Simulated values".
#' @param title title in the convergence plot.
#'
#' @details The function generates the trace plot of simulated values across iterations.
#'
#' @return The plot of simulated values across iterations.
#'
#' @export
conv.plot <- function(data,
                      iter,
                      x.lab = "Iteration number",
                      y.lab = "Simulated values",
                      title = NULL) {
  value <- NULL

  # transform data from wide to long

  wide.dat <- as.data.frame(cbind(data, 0:iter))
  colnames(wide.dat) <- c(colnames(data), "iter")

  # make sure the id column should be factor variable

  long.dat <- melt(wide.dat, id.vars = "iter",
                    measure.vars = colnames(wide.dat)[colnames(wide.dat) != "iter"],
                    variable.name = "var",
                    value.name = "value")

  ggplot(data = long.dat,
         aes(x = as.numeric(iter),
             y = value,
             group = var,
             color = var)) +
    geom_line() +
    xlab(x.lab) +
    ylab(y.lab) +
    ggtitle(title)

}
