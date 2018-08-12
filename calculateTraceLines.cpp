#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List calculateTraceLines(NumericVector theta,
                         NumericVector thr,
                         NumericVector sl,
                         NumericVector wt,
                         int nItems,
                         int nQuad) {
  
  List allTlines(nItems); // list to contain all the trace lines for all items
  
  for(int itemCounter = 0; itemCounter < nItems; ++itemCounter) {
    
    List binomItemList(2); // list of two for the two thresh
    NumericVector pTline(nQuad); // numeric vector containing the trace line of the correct response
    NumericVector qTline(nQuad); // numeric vector containing the trace line of the incorrect response
    
    // calculating tlines
    //for(int qp = 0; qp < nQuad; ++qp) {
      pTline = 1/(1+exp(-sl[itemCounter]*(theta-thr[itemCounter])));
      qTline = 1-pTline;
    //}
    
    double pWeight = wt[itemCounter]; // weight for item x for the correct response
    double qWeight = 0; // 0 for incorect response
    
    List pList(2);
    List qList(2);
    
    pList[0] = pTline;
    pList[1] = pWeight;
    qList[0] = qTline;
    qList[1] = qWeight;
    
    binomItemList[0] = qList;
    binomItemList[1] = pList;
    
    allTlines[itemCounter] = binomItemList;
  }
  return allTlines;
}
