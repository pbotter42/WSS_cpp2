#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]

List calculateRpLik(List allTlines,
           NumericMatrix testeeRp,
           int nItems) {
  
  List rpLikList(testeeRp.nrow());
  
  // start loop for all testees to assign first tline and weight
  for(int testee = 0; testee < testeeRp.nrow(); ++testee) {
    int item1Response = testeeRp(testee, 0); // first item score
    if(item1Response == 0) {
      List testeeList(2); // list containing testee current tLine and wss
      NumericVector testeeTline1 = as<List>(as<List>(as<List>(allTlines[0])[0]))[0]; // tline of the first item for an incorrect response
      std::vector<double> testeeWss1 = as<List>(as<List>(as<List>(allTlines[0])[0]))[1]; // item weight of the first item for an incorrect response
      testeeList[0] = testeeTline1;
      testeeList[1] = testeeWss1;    
      rpLikList[testee] = testeeList;        
    } else {
      List testeeList(2); // list containing testee current tLine and wss
      NumericVector testeeTline1 = as<List>(as<List>(as<List>(allTlines[0])[1]))[0]; // tline of the first item for a correct response
      std::vector<double> testeeWss1 = as<List>(as<List>(as<List>(allTlines[0])[1]))[1]; // weight of the first item for a correct response
      testeeList[0] = testeeTline1;
      testeeList[1] = testeeWss1;    
      rpLikList[testee] = testeeList;
    }
    
    // start items 2:n
    for(int item_n = 1; item_n < testeeRp.ncol(); ++item_n) {
      int item_n_Response = testeeRp(testee, item_n);
      if(item_n_Response == 0) {
        List testeeList(2); // list containing testee current tLine and wss
        
        NumericVector testeePreviousTline = as<List>(as<List>(rpLikList[testee]))[0]; // tLine as of item_n - 1
        std::vector<double> testeePreviousWss = as<List>(as<List>(rpLikList[testee]))[1]; // wss as of item_n - 1
        
        NumericVector currentTline = as<List>(as<List>(as<List>(allTlines[item_n])[0]))[0]; // the nth item tLine for an incorrect response
        std::vector<double> currentWeight = as<List>(as<List>(as<List>(allTlines[item_n])[0]))[1]; // the nth item weight for an incorrect response
        
        NumericVector updatedTline = testeePreviousTline * currentTline;
        std::vector<double> updatedWeight;
        updatedWeight.push_back(testeePreviousWss[0] + currentWeight[0]);
        
        testeeList[0] = updatedTline;
        testeeList[1] = updatedWeight;    
        rpLikList[testee] = testeeList;
        
      } else {
        List testeeList(2); // list containing testee current tLine and wss
        
        NumericVector testeePreviousTline = as<List>(as<List>(rpLikList[testee]))[0]; // tLine as of item_n - 1
        std::vector<double> testeePreviousWss = as<List>(as<List>(rpLikList[testee]))[1]; // wss as of item_n - 1
        
        NumericVector currentTline = as<List>(as<List>(as<List>(allTlines[item_n])[1]))[0]; // the nth item tLine for a correct response
        std::vector<double> currentWeight = as<List>(as<List>(as<List>(allTlines[item_n])[1]))[1]; // the nth item weight for a correct response
        NumericVector updatedTline = testeePreviousTline * currentTline; 
        std::vector<double> updatedWeight;
        updatedWeight.push_back(testeePreviousWss[0] + currentWeight[0]);
        
        testeeList[0] = updatedTline;
        testeeList[1] = updatedWeight;   
        rpLikList[testee] = testeeList;
        
      }
    }
  }
  return rpLikList;
}
