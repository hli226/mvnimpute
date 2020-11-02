#' Function for missing and censored data visualization
#'
#' Draws the plots that graphically presents the percentage of the missing, censored and observed data
#'
#' @param data dataset that contains the missing and censored values
#' @param miss.index missing index indicating the missing data in the original dataset
#' @param censor.index censoring index indicating the censored data in the original dataset
#'
#' @details The function draws the plot that graphically shows the percentage of the missing, censored and observed
#' data in the dataset. Column names of the \code{miss.indx} and \code{censor.indx} matrices should correspond
#' to the respective variable names that are subject to either missing or censoring, and those names should be
#' handled beforehand.
#'
#' @return graph that shows the details of the different type of data in the dataset
#'
#' @export
visual.plot <- function(data,
                        miss.index,
                        censor.index) {

  `%notin%` <- Negate(`%in%`)

  plot.dat <- as.data.frame(data)
  miss.indx <- as.data.frame(miss.index)
  censor.indx <- as.data.frame(censor.index)

  miss.name <- colnames(miss.indx)
  censor.name <- colnames(censor.indx)

  common.name  <- NULL

  for (i in 1:length(miss.name)) {

    if (miss.name[i] %in% censor.name) common.name <- c(common.name, miss.name[i])

  }

  # (1) if all the variable have either missing or censored values
  if (is.null(common.name)) {
    missing.vals <- miss.indx %>%
      gather(key = "key", value = "val") %>%
      mutate(isna = (.data$val == 1)) %>%
      group_by(.data$key) %>%
      mutate(total = n()) %>%
      group_by(.data$key, .data$total, .data$isna) %>%
      summarise(num.isna = n()) %>%
      mutate(pct = .data$num.isna / .data$total * 100)

    censoring.vals <- censor.indx %>%
      gather(key = "key", value = "val") %>%
      mutate(iscensor = (.data$val == 1)) %>%
      group_by(.data$key) %>%
      mutate(total = n()) %>%
      group_by(.data$key, .data$total, .data$iscensor) %>%
      summarise(num.iscensor = n()) %>%
      mutate(pct = .data$num.iscensor / .data$total * 100)

    missing.vals <- data.frame(missing.vals, type = ifelse(missing.vals$isna == TRUE, "0", "1"))
    missing.vals <- missing.vals[, -3]
    colnames(missing.vals) <- c("key", "total", "num", "pct", "type")

    censoring.vals <- censoring.vals[order(censoring.vals[, 1]), ]

    censoring.vals <- data.frame(censoring.vals, type = ifelse(censoring.vals$iscensor == TRUE, "2", "1"))
    censoring.vals <- censoring.vals[, -3]
    colnames(censoring.vals) <- c("key", "total", "num", "pct", "type")

  }

  # if some of the variables have both missing and censored values
  else if (!is.null(common.name)) {

    missing.vals <- miss.indx %>%
      gather(key = "key", value = "val") %>%
      mutate(isna = (.data$val == 1)) %>%
      group_by(.data$key) %>%
      mutate(total = n()) %>%
      group_by(.data$key, .data$total, .data$isna) %>%
      summarise(num.isna = n()) %>%
      mutate(pct = .data$num.isna / .data$total * 100) %>%
      filter(.data$isna == TRUE)

    missing.other <- miss.indx %>%
      gather(key = "key", value = "val") %>%
      mutate(isna = (.data$val == 1)) %>%
      group_by(.data$key) %>%
      mutate(total = n()) %>%
      group_by(.data$key, .data$total, .data$isna) %>%
      summarise(num.isna = n()) %>%
      mutate(miss.pct = .data$num.isna / .data$total * 100) %>%
      filter(.data$isna == FALSE)

    censoring.vals <- censor.indx %>%
      gather(key = "key", value = "val") %>%
      mutate(iscensor = (.data$val == 1)) %>%
      group_by(.data$key) %>%
      mutate(total = n()) %>%
      group_by(.data$key, .data$total, .data$iscensor) %>%
      summarise(num.iscensor = n()) %>%
      mutate(pct = .data$num.iscensor / .data$total * 100) %>%
      filter(!is.na(.data$iscensor))

    missing.vals <- data.frame(missing.vals, type = rep("0", nrow(missing.vals)))
    missing.vals <- missing.vals[, -3]
    colnames(missing.vals) <- c("key", "total", "num", "pct", "type")

    missing.other <- data.frame(missing.other, type = rep("1", nrow(missing.other)))
    missing.other <- missing.other[, -3]
    colnames(missing.other) <- c("key", "total", "num", "pct", "type")

    censoring.vals <- data.frame(censoring.vals, type = ifelse(censoring.vals$iscensor == TRUE, "2", "1"))
    censoring.vals <- censoring.vals[, -3]
    colnames(censoring.vals) <- c("key", "total", "num", "pct", "type")

    missing.other <- missing.other %>%
      filter(missing.other$key %notin% censoring.vals$key)

    missing.vals <- rbind(missing.vals, missing.other)

    censoring.vals <- censoring.vals[order(censoring.vals[, 1]), ]


  }

  no.miss.censor <- colnames(data)[colnames(data) %notin% c(unique(missing.vals$key), unique(censoring.vals$key))]

  if (length(no.miss.censor) == 0) {
    complete.dat <- rbind(missing.vals,
                          censoring.vals)
  } else {

    no.miss.censor.mat <- NULL

    for (i in 1:nrow(as.data.frame(no.miss.censor))) { # NOTE: as.data.frame transform vector to a column matrix

      no.miss.censor.mat <- rbind(no.miss.censor.mat,
                                  data.frame(key = no.miss.censor[i],
                                             total = nrow(data),
                                             num = nrow(data),
                                             pct = 100,
                                             type = 1))

    }

    colnames(no.miss.censor.mat) <- c("key", "total", "num", "pct", "type")

    complete.dat <- rbind(missing.vals,
                          censoring.vals,
                          no.miss.censor.mat)
  }

  complete.dat <- complete.dat %>%
    arrange(.data$key, .data$type)

  complete.dat %>%
    ggplot() +
    geom_bar(aes(x = .data$key,
                 y = as.numeric(.data$pct),
                 fill = .data$type), stat = "identity") +
    scale_fill_manual(name = "",
                      values = c("tomato3", "steelblue", "lightgreen"),
                      labels = c("Missing", "Observed", "Censored")) +
    coord_flip() +
    xlab("Variable") +
    ylab("Percentage") +
    ggtitle("Percentage of missing/censored values in each variable")

  # return(complete.dat)

}
