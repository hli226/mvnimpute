#' Marginal density plots function
#'
#' Draws the marginal density plots for each variable
#'
#' @param data dataset containing all the variables
#' @param title title of the generated plots
#'
#' @return Marginal density plot for each variable in the dataset
#'
#' @export
marg.plot <- function(data, title = NULL) {

  # draw marginal plots to visual distribution of each marginal
  if (!is.null(ncol(data))) {

    for (i in 1:ncol(data)) {

      plot(density(na.omit(data[, i])), main = title[i])

    }

  } else {

    plot(density(na.omit(data)), main = title)

  }

}
