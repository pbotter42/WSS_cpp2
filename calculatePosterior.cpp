#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calculatePosterior(NumericVector dist, NumericMatrix summedTraceLines) {
  
  NumericMatrix out(summedTraceLines.nrow(), summedTraceLines.ncol());
  for(int i = 0; i < summedTraceLines.nrow(); ++i) {
  out.row(i) = dist * summedTraceLines.row(i);
  }
  
  return out;
}
