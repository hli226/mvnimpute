
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mvnimpute

<!-- badges: start -->

<!-- badges: end -->

The goal of mvnimpute is to implement multiple imputation to the data
when there are both missing and censored values (a single variable can
have missing and censoring simultaneously; or it can have either one of
missing or censored values)

## Installation

You can the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("yuebanfengqing/mvnimpute")
```

## Installation

You can install the version of mvnimpute in development from
[GitHub](https://github.com) with:

``` r
# install.packages("devtools")
devtools::install_packages_github("yuebanfengqing/mvnimpute")
```

add `built_vignettes = TRUE` argument to include vignette in the
downloaded package. You have to install the development package
**devtools** for installing packages from GitHub.

## Basic functions

It has 10 functions including

`data.gen`: generates multivariate normal data with missing and
censoring information.

`calcu.param`: calculates the complete-case and available-case
parameters.

`visual.plot`: draw plot showing percentage of missing, censored, and
observed values.

`marg.plot`: draw density plot for each marginal variable.

`initial.impute`: implements single imputation approach to makeup
incomplete data.

`multiple.impute`: implements multiple imputation for missing and
censored data, have to perform on the complete data.

`conv.plot`: draws trace plots of the parameters from the multiple
imputation.

`avg.plot`: draws trace plots of the averaged values of the parameters
from the multiple imputation.

`calcu.acf`: calculates the autocorrelation values and draws acf plots.

`nhanes.dat`: 2011 - 2016 NHANES demographics, PCB measurements, and
heavy metal pollutants data.

For detailed instructions on using the package and a walk-through
example, please refer to the vignette downloaded with the package.
