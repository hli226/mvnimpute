#' Function for missing and censored data visualization
#'
#' Draws the plots that graphically presents the percentage of the missing, censored and observed data
#'
#' @param data.ind dataset with missing and censored values
#' @param title a
#'
#' @details The function draws the plot that graphically shows the percentage of the missing, censored and observed
#' data in the dataset. Column names of the \code{miss.indx} and \code{censor.indx} matrices should correspond
#' to the respective variable names that are subject to either missing or censoring, and those names should be
#' handled beforehand.
#'
#' @return The plot that shows the details of the different type of data in the dataset
#'
#' @export
visualization.plot <- function(
  data.ind,
  title = NULL) {

  n <- nrow(data.ind)
  p <- ncol(data.ind)

  data.t <- matrix(NA, nrow = p, ncol = 3)
  colnames(data.t) <- c("Observed", "Missing", "Censored")

  for (i in 1:p) {

    summ <- table(data.ind[, i])
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

  new.data <- data.frame(id = paste0("y", 1:p), data.t)
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
