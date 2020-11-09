#' Function for missing and censored data visualization
#'
#' Draws the plots that graphically presents the percentage of the missing, censored and observed data
#'
#' @param data dataset that contains the missing and censored values
#' @param miss.index missing index indicating the missing data in the original dataset
#' @param censor.index censoring index indicating the censored data in the original dataset
#' @param title the title of the plot. Default is NULL.
#'
#' @details The function draws the plot that graphically shows the percentage of the missing, censored and observed
#' data in the dataset. Column names of the \code{miss.indx} and \code{censor.indx} matrices should correspond
#' to the respective variable names that are subject to either missing or censoring, and those names should be
#' handled beforehand.
#'
#' @return graph that shows the details of the different type of data in the dataset
#'
#' @export
visual.plot <- function(
  data,
  miss.index = NULL,
  censor.index = NULL,
  title = NULL) {

  `%notin%` <- Negate(`%in%`)

  # 1. data only has missing values
  if (!is.null(miss.index) & is.null(censor.index)) {
    # either missing or censoring index is null matrix
    plot.dat <- data
    miss.indx <- as.data.frame(miss.index)

    missing.vals <- miss.indx %>%
      gather(key = "key", value = "val") %>%
      mutate(isna = (.data$val == 1)) %>%
      group_by(.data$key) %>%
      mutate(total = n()) %>%
      group_by(.data$key, .data$total, .data$isna) %>%
      summarise(num.isna = n()) %>%
      mutate(pct = .data$num.isna / .data$total * 100)

    no.missing <- colnames(plot.dat)[colnames(plot.dat) %notin% unique(missing.vals$key)]

    if (length(no.missing) > 0) {

      no.missingfalse <- data.frame(
        key = no.missing,
        total = nrow(plot.dat),
        isna = FALSE,
        num.isna = nrow(plot.dat),
        pct = 100
      )

      no.missingtrue <- data.frame(
        key = no.missing,
        total = nrow(plot.dat),
        isna = TRUE,
        num.isna = 0,
        pct = 0
      )

      missing.vals <- rbind(missing.vals, no.missingfalse, no.missingtrue)

    }

    levels <- (missing.vals %>% filter(isna == TRUE) %>% arrange(desc(pct)))$key

    missing.vals %>%
      ggplot() +
      geom_bar(aes(x = reorder(.data$key, desc(.data$pct)),
                   y = .data$pct, fill = .data$isna),
               stat = "identity") +
      scale_x_discrete(limits = levels,
                       labels = paste(
                         levels,
                         " (",
                         round((missing.vals %>% filter(isna == TRUE) %>% arrange(desc(pct)))$pct, 2),
                         "%)", sep = ""
                       )) +
      scale_fill_manual(name = "",
                        values = c("steelblue", "tomato3"),
                        labels = c("Observed", "Missing")) +
      coord_flip() +
      labs(
        x = "Variable",
        y = "Percentage") +
      ggtitle(title)

  }

  # 2. data only has censored values
  else if (is.null(miss.index) & !is.null(censor.index)) {

    plot.dat <- data
    censor.indx <- as.data.frame(censor.index)

    censoring.vals <- censor.indx %>%
      gather(key = "key", value = "val") %>%
      mutate(iscensor = (.data$val == 1)) %>%
      group_by(.data$key) %>%
      mutate(total = n()) %>%
      group_by(.data$key, .data$total, .data$iscensor) %>%
      summarise(num.iscensor = n()) %>%
      mutate(pct = .data$num.iscensor / .data$total * 100)

    no.censoring <- colnames(plot.dat)[colnames(plot.dat) %notin% unique(censoring.vals$key)]

    if (length(no.censoring) > 0) {

      no.censoringfalse <- data.frame(
        key = no.censoring,
        total = nrow(plot.dat),
        iscensor = FALSE,
        num.iscensor = nrow(plot.dat),
        pct = 100
      )

      no.censoringtrue <- data.frame(
        key = no.censoring,
        total = nrow(plot.dat),
        iscensor = TRUE,
        num.iscensor = 0,
        pct = 0
      )

      censoring.vals <- rbind(censoring.vals, no.censoringfalse, no.censoringtrue)

    }

    levels <- (censoring.vals %>% filter(iscensor == TRUE) %>% arrange(desc(pct)))$key

    censoring.vals %>%
      ggplot() +
      geom_bar(aes(x = reorder(.data$key, desc(.data$pct)),
                   y = .data$pct, fill = .data$iscensor),
               stat = "identity") +
      scale_x_discrete(limits = levels,
                       labels = paste(
                         levels,
                         " (",
                         round((censoring.vals %>% filter(iscensor == TRUE) %>% arrange(desc(pct)))$pct, 2),
                         "%)", sep = ""
                       )) +
      scale_fill_manual(name = "",
                        values = c("steelblue", "lightgreen"),
                        labels = c("Observed", "Censored")) +
      coord_flip() +
      labs(
        x = "Variable",
        y = "Percentage") +
      ggtitle(title)

  }

  else if (!is.null(miss.index) & !is.null(censor.index)) {

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

      observing <- colnames(plot.dat)[colnames(plot.dat) %notin% c(unique(missing.vals$key), unique(censoring.vals$key))]

      if (length(observing) > 0) {

      observing.vals <- data.frame(
        key = observing,
        total = nrow(plot.dat),
        num = rep(nrow(plot.dat), length(observing)),
        pct = 100,
        type = 1
      )
      }
      complete.dat <- rbind(missing.vals, censoring.vals, observing.vals)

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
        mutate(pct = .data$num.isna / .data$total * 100) %>%
        filter(.data$isna == FALSE)

      no.missing <- missing.other$key[missing.other$key %notin% missing.vals$key]

      if (length(no.missing) > 0) {

        missing.0 <- data.frame(
          key = no.missing,
          total = nrow(plot.dat),
          isna = rep(TRUE, length(no.missing)),
          num.isna = rep(0, length(no.missing)),
          pct = rep(0, length(no.missing))
        )

        missing.vals <- rbind(missing.vals, missing.0)
      }

      censoring.vals <- censor.indx %>%
        gather(key = "key", value = "val") %>%
        mutate(iscensor = (.data$val == 1)) %>%
        group_by(.data$key) %>%
        mutate(total = n()) %>%
        group_by(.data$key, .data$total, .data$iscensor) %>%
        summarise(num.iscensor = n()) %>%
        mutate(pct = .data$num.iscensor / .data$total * 100) %>%
        filter(.data$iscensor == TRUE)

      censoring.other <- censor.indx %>%
        gather(key = "key", value = "val") %>%
        mutate(iscensor = (.data$val == 1)) %>%
        group_by(.data$key) %>%
        mutate(total = n()) %>%
        group_by(.data$key, .data$total, .data$iscensor) %>%
        summarise(num.iscensor = n()) %>%
        mutate(pct = .data$num.iscensor / .data$total * 100) %>%
        filter(.data$iscensor == FALSE)

      no.censoring <- censoring.other$key[censoring.other$key %notin% censoring.vals$key]

      if (length(no.censoring) > 0) {

        censoring.0 <- data.frame(
          key = no.censoring,
          total = nrow(plot.dat),
          iscensor = rep(TRUE, length(no.censoring)),
          num.iscensor = rep(0, length(no.censoring)),
          pct = rep(0, length(no.censoring))
        )

        censoring.vals <- rbind(censoring.vals, censoring.0)

      }

      missing.vals <- data.frame(missing.vals, type = rep("0", nrow(missing.vals)))
      missing.vals <- missing.vals[, -3]
      colnames(missing.vals) <- c("key", "total", "num", "pct", "type")
      missing.vals <- missing.vals[order(missing.vals$key), ]

      censoring.vals <- data.frame(censoring.vals, type = ifelse(censoring.vals$iscensor == TRUE, "2", "1"))
      censoring.vals <- censoring.vals[, -3]
      colnames(censoring.vals) <- c("key", "total", "num", "pct", "type")
      censoring.vals <- censoring.vals[order(censoring.vals$key), ]

      observing.vals <- data.frame(
        key = missing.vals$key,
        total = rep(nrow(plot.dat), length(missing.vals$key)),
        num = nrow(plot.dat) - missing.vals$num - censoring.vals$num,
        pct = (nrow(plot.dat) - missing.vals$num - censoring.vals$num) / nrow(plot.dat),
        type = rep(1, length(missing.vals$key))
      )

      complete.dat <- rbind(missing.vals,
                            censoring.vals,
                            observing.vals)

      complete.dat$pct <- complete.dat$num / complete.dat$total * 100

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
      ggtitle(title)

  }

}
