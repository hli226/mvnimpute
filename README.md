
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mvnimpute

<!-- badges: start -->
<!-- badges: end -->

The goal of **mvnimpute** package is to implement multiple imputation to
the data when there are both missing and censored values (a single
variable can have missing and censoring simultaneously; or it can have
either only missing or censored values). You can either install this
package directly from GitHub (R &gt;= 3.6.0) or from local source file.

## Installation

### From GitHub

**NOTE: To correctly install the package from GitHub, you have to update
the R software to at least version 3.6.0**. You can install
**mvnimpute** package in development version from
[GitHub](https://github.com) with:

``` r
# install.packages("devtools")
devtools::install_github("yuebanfengqing/mvnimpute")
```

add `build_vignettes = TRUE` argument to include vignette in the
downloaded package. You have to install the development package
**devtools** for installing packages from GitHub. The packages that
**mvnimpute** depends on will be automatically downloaded and installed.

### From local file

For R version number less than **3.5.0**, you have to install the
corresponding Rtools from
<https://cran.r-project.org/bin/windows/Rtools/history.html> to the
system PATH before running the following code.

Please email to <hli226@uic.edu> for the latest compiled package source
file. For installing package from local source file, you have to
manually install the dependencies first as

``` r
## install dependencies
install.packages("ggplot2")
install.packages("truncnorm")
install.packages("mvtnorm")
install.packages("reshape2")
install.packages("LaplacesDemon")
install.packages("dplyr")
install.packages("magrittr")
install.packages("tidyr")
install.packages("rlang")
```

Then the package can be installed with replacing the *path-of-file* by
your local path that stores the package source file
(**mvnimpute\_0.0.0.9000.tar.gz**) as

``` r
## install mvnimpute package from local source file
install.packages("path-of-file/mvnimpute_0.0.0.9000.tar.gz", repos = NULL, type = "source")
```

## Basic functions

It has 10 functions including

`data.gen`: generates multivariate normal data with missing and
censoring information.

`calcu.param`: calculates the complete-case and available-case
parameters.

`visual.plot`: draws plot showing percentage of missing, censored, and
observed values.

`marg.plot`: draws density plot for each marginal variable.

`initial.impute`: implements single imputation approach to makeup
incomplete data.

`multiple.impute`: implements multiple imputation for missing and
censored data, have to perform on the complete data.

`conv.plot`: draws trace plots of the parameters from the multiple
imputation.

`avg.plot`: draws trace plots of the averaged values of the parameters
from the multiple imputation.

`calcu.acf`: calculates the autocorrelation values and draws ACF plots.

`nhanes.dat`: 2011 - 2016 NHANES demographics, PCB measurements, and
heavy metal pollutants data.

For detailed instructions on using the package and a walk-through
example, please refer to the vignette downloaded with the package. Once
the package is installed, the vignette can be accessed using the
function `browseVignettes("mvnimpute")`, and it will open in an external
viewer.
