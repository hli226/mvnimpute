#' Convergence plot function
#'
#' This function generates the convergence plots for the parameter values
#'
#' @param dat The dataset containing the values to be plotted
#' @param iter Iteration number of the simulation
#' @param x.lab Label of the x axis in the generated plot, default is set to be "Iteration number"
#' @param y.lab Label of the y axis in the generated plot, default is set to be "Simulated values"
#' @param title Title in the convergence plot
#'
#' @return Convergence plot of your desired variable
#'
#' @export
conv.plot <- function(dat,
                      iter,
                      x.lab = "Iteration Number",
                      y.lab = "Simulated values",
                      title = NULL) {
  value <- NULL

  # transform data from wide to long

  wide.dat <- as.data.frame(cbind(dat, 0:iter))
  colnames(wide.dat) <- c(colnames(dat), "iter")

  # make sure the id column should be factor variable
  # wide.mean$iter <- as.factor(wide.mean$iter)

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
