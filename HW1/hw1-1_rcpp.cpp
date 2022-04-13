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

double g(const double x){
  double value;
  value = std::log(x)/(1+x);
  return value;
}

double gprime1(const double x){
  double value;
  value = (1+(1/x)-std::log(x))/(std::pow((1+x),2));
  return value;
}

double gprime2(const double x){
  double value;
  value = (-1/(std::pow(x,2)+std::pow(x,3)))-2*(1+(1/x)-std::log(x))/(std::pow((1+x),3));
  return value;
}

// [[Rcpp::export]]
double newton(double x, const double epsilon){
  int i = 0;
  double diff = 1.0;
  
  while(std::fabs(diff) > epsilon){
    diff = -gprime1(x)/gprime2(x);
    x = x + diff;
    i = i + 1;
    Rprintf("%.2d %.3f\n",i,x);
  }
  
  return x;
}
