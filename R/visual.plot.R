#' Function for missing and censored data visualization
#'
#' Draws plot that graphically show the percentages of the missing, censored and observed data
#'
#' @param data.indicator matrix that contains the data type indicators of the original data
#' @param title title of the generated plot, default is set to "Summary plot"
#'
#' @details The function draws the plot that graphically shows the percentages of the missing, censored and observed
#' data in the dataset. \code{data.indicator} should be a matrix containing the data type indicators as generated in the
#' data preparation step. 0 for missing values, 1 for obsrved values, and 2 for censored values. \code{title} is the title
#' of the generated plot.
#'
#' @return The plot that shows the details of the different type of data in the dataset
#'
#' @export
visual.plot <- function(
  data.indicator,
  title = "Summary plot") {

  c.names <- colnames(data.indicator) # get column names of the data
  n <- nrow(data.indicator)
  p <- ncol(data.indicator)

  data.t <- matrix(NA, nrow = p, ncol = 3)
  colnames(data.t) <- c("Observed", "Missing", "Censored")

  for (i in 1:p) {

    summ <- table(data.indicator[, i])
    names.summ <- names(summ)

    if ("0" %in% names.summ) {
      data.t[i, "Missing"] <- summ[names.summ == 0]
    } else {
      data.t[i, "Missing"] <- 0
    }

    if ("1" %in% names.summ) {
      data.t[i, "Observed"] <- summ[names.summ == 1]
    } else {
      data.t[i, "Observed"] <- 0
    }

    if ("2" %in% names.summ) {
      data.t[i, "Censored"] <- summ[names.summ == 2]
    } else {
      data.t[i, "Censored"] <- 0
    }

  }

  new.data <- data.frame(id = c.names, data.t)
  new.data$Observed <- new.data$Observed/n * 100
  new.data$Missing <- new.data$Missing/n * 100
  new.data$Censored <- new.data$Censored/n * 100

  percent.mat <- melt(new.data, id.vars = "id")

  ggplot(percent.mat) +
    geom_bar(aes(x = .data$id, y = .data$value, fill = .data$variable),
             stat = "identity") +
    coord_flip() +
    labs(x = "Variable", y = "Percentage (%)") +
    ggtitle(title)
}
