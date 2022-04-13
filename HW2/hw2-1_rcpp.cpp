#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

double g(NumericVector data, const double x){
  double value = 0.0;
  int i;
  int n = data.size();
  
  value = -n * log(3.141592);
  for(i = 0; i < n; i++){
    value += -log(1 + std::pow((data[i]-x),2.0));
  }
  
  return value;
}

double gprime1(NumericVector data, const double x){
  double value = 0.0;
  int i;
  int n = data.size();
  
  for(i = 0; i < n; i++){
    value += 2.0 * (data[i]-x)/(1.0 + std::pow((data[i]-x),2.0));
  }
  
  return value;
}

double gprime2(NumericVector data, const double x){
  double value = 0.0;
  int i;
  int n = data.size();
  
  for(i = 0; i < n; i++){
    value += 2.0 * (std::pow((data[i]-x),2.0)-2.0)/(std::pow((1.0+std::pow((data[i]-x),2.0)),2.0));
  }
  
  return value;
}

// [[Rcpp::export]]
double newton_cauchy(NumericVector data, double x, const double epsilon){
  int i = 0;
  double diff = 1.0;
  double value;
  
  while(std::fabs(diff) > epsilon){
    diff = -gprime1(data,x)/gprime2(data,x);
    x = x + diff;
    i = i + 1;
    value = g(data,x);
    Rprintf("%.2d %.3f %.3f\n",i,x,value);
  }
  
  return x;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//