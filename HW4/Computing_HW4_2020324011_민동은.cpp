#include <Rcpp.h>
using namespace Rcpp;

// Question 5.3 (a)
double f(double x, double mu) {
  double value;
  value = sqrt(7/(18*M_PI)) * exp(-7*std::pow((x-mu),2)/18) / (2*M_PI*(1+std::pow((mu-5)/2, 2)));
  return value;
}

// [[Rcpp::export]]
double trapezoid(double x, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sumf = 0.0;
  for (i=1; i<n; i++){
    sumf += f(x,a+i*h);
    }
  return h*f(x,a)/2 + h*sumf + h*f(x,b)/2;
}

// Question 5.3 (b)
// [[Rcpp::export]]
double riemann(double x, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sumf = 0.0;
  for (i=0; i<n; i++){
    sumf += f(x,a+i*h);
  }
  return h*sumf;
}

// [[Rcpp::export]]
// n is an even number!!
double simpson(double x, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sumf = 0.0;
  for (i=1; i<(n/2)+1; i++){
    sumf += f(x, a+(2*i-2)*h) + 4*f(x, a+(2*i-1)*h) + f(x, a+(2*i*h));
  }
  return (h/3)*sumf;
}


// Question 5.3 (c)
double transf(double x, double u){
  double value = f(x, log(u/(1-u))) / u / (1 - u);
  return value;
}

// [[Rcpp::export]]
double criemann(double x, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sumf = 0.0;
  for (i=0; i<n; i++){
    sumf += transf(x,a+i*h);
  }
  return h*sumf;
}

// [[Rcpp::export]]
double csimpson(double x, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sumf = 0.0;
  for (i=1; i<(n/2)+1; i++){
    sumf += transf(x, a+(2*i-2)*h) + 4*transf(x, a+(2*i-1)*h) + transf(x, a+(2*i*h));
  }
  return (h/3)*sumf;
}


// Question 5.3 (d)
double transf2(double x, double u){
  double value = f(x, (1/u)) * std::pow((1/u), 2);
  return value;
}

// [[Rcpp::export]]
double dsimpson(double x, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sumf = 0.0;
  for (i=1; i<(n/2)+1; i++){
    sumf += transf2(x, a+(2*i-2)*h) + 4*transf2(x, a+(2*i-1)*h) + transf2(x, a+(2*i*h));
  }
  return (h/3)*sumf;
}


// Question 5.4
// [[Rcpp::export]]
double g(double x) {
  double value = 1/x;
  return value;
}

// [[Rcpp::export]]
double gtrapezoid(double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sumg = 0.0;
  for (i=1; i<n; i++){
    sumg += g(a+i*h);
  }
  return h*g(a)/2 + h*sumg + h*g(b)/2;
}

