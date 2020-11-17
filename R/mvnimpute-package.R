#' mvnimpute: Implementation of multiple imputation involving the missing and censored values
#'
#' The mvnimpute package implements multiple imputation for concurrently handling missing and censored values based on the joint normal
#' model assumption.
#'
#' @docType package
#'
#' @name mvnimpute-package
#'
#' @author
#'  Hesen Li
NULL
# various imports
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats quantile
#' @importFrom stats rexp
#' @importFrom stats na.omit
#' @importFrom stats var
#' @importFrom stats cor
#' @importFrom stats cov
#' @importFrom stats rnorm
#' @importFrom stats complete.cases
#' @importFrom stats density
#' @importFrom stats runif
#' @importFrom stats reorder
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom LaplacesDemon rinvwishart
#' @importFrom truncnorm rtruncnorm
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom reshape2 melt
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr n
#' @importFrom plyr desc
#' @importFrom rlang .data
NULL

