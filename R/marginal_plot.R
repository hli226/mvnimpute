#' Marginal plots function
#'
#' Draws the marginal density plots for each variable
#'
#' @param dat Dataset containing all the variables
#' @param title Titles of the generated plots
#'
#' @return Marginal density plots for each variable in the dataset
#'
#' @export
marg.plot <- function(dat, title = NULL) {

  # draw marginal plots to visual distribution of each marginal
  if (!is.null(ncol(dat))) {

    for (i in 1:ncol(dat)) {

      plot(density(na.omit(dat[, i])), main = title[i])

    }

  } else {

    plot(density(na.omit(dat)), main = title)

  }

}
