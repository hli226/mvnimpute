#' Marginal density plots function
#'
#' This function draws the marginal density plots for all variables.
#'
#' @param data.mat data matrix including all the variables.
#' @param title title of the generated plots.
#'
#' @return Marginal density plot for each variable in the dataset.
#'
#' @export
marg.plot <- function(data.mat, title = NULL) {

  # draw marginal plots to visual distribution of each marginal
  if (!is.null(ncol(data.mat))) {

    for (i in 1:ncol(data.mat)) {

      plot(density(na.omit(data.mat[, i])), main = title[i])

    }

  } else {

    plot(density(na.omit(data.mat)), main = title)

  }

}
