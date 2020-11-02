---
title: mvnimpute package NEWS
---
# mvnimpute 0.0.0.9000

* Create the package for implementing multiple imputation.
* Currently, the *dat.gen* function only supports generation of the data with certain types of missing and censored values
* Currently this package deals with the missing and censored data under multivariate normality assumption.
* For single imputation approach applied in the *initial.impute* function, it currently only supports the stochastic regression imputation, need to extend this function to incorporate more SI approaches in the future.
* For multiple imputation in the *multiple.impute* function, it uses the Normal-Inverse-Wishart conjugate prior distributions for performing the algorithm. Next will extend to more priors.
