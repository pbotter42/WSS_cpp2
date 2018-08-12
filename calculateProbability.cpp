#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int timesTwo(NumericVector x) {
  return x * 2;
}



