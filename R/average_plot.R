#' Average values plot function
#'
#' This function generates average simulated values of the parameter values
#'
#' @param dat Dataset containing the simulated values
#' @param iter Iteration number of the simulation
#' @param x.lab Label of the x axis in the generated plot, default is set to be "Iteration number"
#' @param y.lab Label of the y axis in the generated plot, default is set to be "Average of simulated values"
#' @param title Title of the generated plot
#' @param details Logical variable to specify whether the average simulate values are returned, default is set to TRUE
#'
#' @export
avg.plot <- function(dat,
                     iter,
                     x.lab = "Iteration number",
                     y.lab = "Average of simulated value",
                     title = NULL,
                     details = TRUE) {

  # compute the cumulative averages of the parameter values to the current simulation
  avg.param <- matrix(nrow = iter + 1, ncol = ncol(dat))
  avg.param[1, ] <- dat[1, ]

  for (i in 2:nrow(avg.param)) {

    avg.param[i, ] <- apply(dat[1:i, ], 2, mean)

  }

  # transform data from wide to long

  wide.dat <- as.data.frame(cbind(avg.param, 0:iter))
  colnames(wide.dat) <- c(colnames(dat), "iter")

  # make sure the id column should be factor variable
  # wide.mean$iter <- as.factor(wide.mean$iter)

  long.dat <- melt(wide.dat, id.vars = "iter",
                   measure.vars = colnames(wide.dat)[colnames(wide.dat) != "iter"],
                   variable.name = "var",
                   value.name = "value")

  print(ggplot(data = long.dat,
         aes(x = as.numeric(iter),
             y = value,
             group = var,
             color = var)) +
    geom_line() +
    xlab(x.lab) +
    ylab(y.lab) +
    ggtitle(title))

  return(long.dat)

}
