
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mvnimpute

<!-- badges: start -->

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fhli226%2Fmvnimpute&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=visits&edge_flat=false)](https://hits.seeyoufarm.com)
![Visitor
count](https://shields-io-visitor-counter.herokuapp.com/badge?page=https://github.com/hli226/mvnimpute&color=pink&label=visits)
![GitHub release (latest by
date)](https://img.shields.io/github/v/release/hli226/mvnimpute)
![GitHub all
releases](https://img.shields.io/github/downloads/hli226/mvnimpute/total?color=green)
![GitHub release (by
tag)](https://img.shields.io/github/downloads/hli226/mvnimpute/latest/total)
![GitHub top
language](https://img.shields.io/github/languages/top/hli226/mvnimpute?color=red)
![GitHub language
count](https://img.shields.io/github/languages/count/hli226/mvnimpute?color=purple)
![GitHub R package
version](https://img.shields.io/github/r-package/v/hli226/mvnimpute)

<!-- badges: end -->

The goal of **mvnimpute** package is to implement multiple imputation to
the data when there are both missing and censored values (a single
variable can have both missing and censored values simultaneously; or it
can have either only missing or censored values). An example of
application of this package is for imputing the NHANES laboratory
measurement data that are subject to limits of detection (LODs).

## Installation

For Windows users, the Rtools for building R packages have to be
installed according to your R version from
<https://cran.r-project.org/bin/windows/Rtools/history.html>.

### From GitHub

**NOTE: Some of the packages this package depends on may require the
latest version of R, it is recommended to update your R software to the
latest version**. You can install **mvnimpute** package in its
development version from [GitHub](https://github.com) with:

#### For first-time users

``` r
# install the development package devtools for installing packages from GitHub
install.packages("devtools")

# install mvnimpute package with vignette from GitHub
devtools::install_github("hli226/mvnimpute", build_vignettes = TRUE)
```

`build_vignettes = TRUE` argument is added for including the vignettes,
which give the step-by-step instructions on how to use this package
using an artificial example and a real data example. You have to install
the development package **devtools** for installing packages from
GitHub. The packages that **mvnimpute** depends on will be automatically
downloaded and installed.

## Basic functions

It has 8 functions including

`data.generation`: generates multivariate normal data with missing and
censored values.

`param.calc`: calculates the complete-case and available-case
parameters.

`visual.plot`: draws plot showing percentages of missing, censored, and
observed values.

`marg.plot`: draws density plot for each marginal variable.

`multiple.imputation`: implements multiple imputation for missing and
censored data, have to perform on the complete data.

`conv.plot`: draws trace plot of the parameters from the multiple
imputation.

`avg.plot`: draws trace plot of the averaged values of the parameters
from the multiple imputation.

`acf.calc`: calculates the autocorrelation values and draws ACF plots.

`nhanes.dat`: 2011 - 2016 NHANES demographics, PCB measurements, and
heavy metal pollutants data.

For detailed instructions on using the package and a walk-through
example, please refer to the vignette downloaded with the package. Once
the package is installed, the vignette can be accessed using the
function `browseVignettes("mvnimpute")`, and it will open in an external
web browser.

## Acknowlegements

This package is based on the work supported by the National Institute of
Environmental Health Sciences under grant 1R01ES028790-01.
