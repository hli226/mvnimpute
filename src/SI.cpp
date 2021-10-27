// [[Rcpp::depends(RcppDist)]]
#include <RcppDist.h>
using namespace Rcpp;

//' Calculate CC parameters
//'
//' @param data a list
//' @export
// [[Rcpp::export]]
List param_calc(const List data) {

  NumericMatrix lval = data[0];
  NumericMatrix rval = data[1];

  const int n = lval.nrow(), p = lval.ncol();

  // NumericMatrix obs_indx(n, p);

  // CC imputation
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < p; j++) {
  //     if (lval(i, j) == rval(i, j))  // if the values is observed
  //       obs_indx(i, j) = 1;
  //     else if (lval(i, j) != rval(i, j) && lval(i, j) == -10000 && rval(i, j) < 10000) // left censored
  //       obs_indx(i, j) = 1;
  //     else if (lval(i, j) != rval(i, j) && lval(i, j) > -10000 && rval(i, j) == 10000) // right censored
  //       obs_indx(i, j) = 1;
  //   }
  // }
  //
  // // count how many CC cases
  // NumericVector count_cc(n);
  //
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < p; j++) {
  //     if (obs_indx(i, j) == 1)
  //       count_cc(i)++;
  //   }
  // }
  //
  // int cc_case = 0;
  // for (int i = 0; i < n; i++) {
  //   if (count_cc(i) == p)
  //     cc_case++;
  // }
  // NumericMatrix CC_dat(cc_case, p);
  //
  // for (int i = 0, count = 0; i < n; i++) {
  //   if (count_cc(i) == p) {
  //     for (int j = 0; j < p; j++) {
  //       if (lval(i, j) == -10000)
  //         CC_dat(count, j) = rval(i, j);   // left censored
  //       else if (rval(i, j) == 10000)
  //         CC_dat(count, j) = lval(i, j);   //  righe censored
  //     }
  //     count++;
  //   }
  // }
  // // calculate CC mean and variance
  // NumericVector CC_mean(p), CC_var(p);
  //
  // for (int j = 0; j < p; j++) {
  //   CC_mean(j) = mean(CC_dat(_, j));
  //   CC_var(j) = var(CC_dat(_, j));
  // }

  NumericVector not_miss(p);  // count number of not missing for each variable
  for (int j = 0; j < p; j++) {
    for (int i = 0; i < n; i++) {
      if (lval(i, j) != rval(i, j) && abs(lval(i, j)) == abs(rval(i, j)))
        continue;
      else
        not_miss(j)++;
    }
  }

  NumericVector AC_mean(p);
  NumericVector AC_var(p);
  // AC imputation
  for (int j = 0; j < p; j++) {
    int len_var = not_miss(j);
    NumericVector obs_val(len_var);

    for (int i = 0, count = 0; i < n; i++) {
      if (lval(i, j) != rval(i, j) && abs(lval(i, j)) == abs(rval(i, j)))
        continue;
      else {
        if (lval(i, j) > -10000 && rval(i, j) == 10000)
          obs_val(count) = lval(i, j);
        else if (lval(i, j) == -10000 && rval(i, j) < 10000)
          obs_val(count) = rval(i, j);
        else
          obs_val(count) = lval(i, j);
      }
    }
    AC_mean(j) = mean(obs_val);
    AC_var(j) = var(obs_val);
  }

  return List::create(
    // Named("CC.mean") = CC_mean,
    // Named("CC.var") = CC_var
    // not_miss

    Named("CC.mean") = AC_mean,
    Named("CC.var") = AC_var
  );

}

//' Single imputation
//'
//' @param data a list
//' @export
// [[Rcpp::export]]
NumericMatrix single_imputation(List data) {

  NumericVector CC_mean = param_calc(data)["CC.mean"];
  NumericVector CC_var = param_calc(data)["CC.var"];

  NumericMatrix lval = data[0], rval = data[1];

  // number of observations and variables
  const int n = lval.nrow(), p = lval.ncol();

  // SI data
  NumericMatrix SI_dat(n, p);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      // (1) observed values
      if (lval(i, j) == rval(i, j))
        SI_dat(i, j) = lval(i, j);
      // (2) missing and censored values
      else
        SI_dat(i, j) = R::rnorm(CC_mean(j), sqrt(CC_var(j)));
      // (3) censored values
      // else {
      //   SI_dat(i, j) = r_truncnorm(CC_mean(j), sqrt(CC_var(j)), lval(i, j), rval(i, j));
      // }
    }
  }

  // return SI_dat;
  return SI_dat;

}
