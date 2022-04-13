// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// question 1-a
// [[Rcpp::export]]
double g(arma::dvec x, double theta){
  double value;
  value = arma::sum(-log(M_PI)-log(1+arma::pow((x-theta),2)));
  return value;
}

double gprime1(arma::dvec x, double theta){
  double value;
  value = arma::sum(2*(x-theta)/(1+arma::pow((x-theta),2)));
  return value;
}

double gprime2(arma::dvec x, double theta){
  double value;
  value = -arma::sum(2*(1-arma::pow((x-theta),2))/(1+arma::pow((x-theta),2)));
  return value;
}

// [[Rcpp::export]]
double newton(arma::dvec x, double init, double eps){
  int i = 0;
  double theta = init;
  while(std::fabs(gprime1(x, theta)) > eps){
    theta = theta - gprime1(x, theta)/gprime2(x, theta);
    i = i + 1;
    // Rprintf("%.2d %.3f\n",i,theta);
  }
  return theta;
}


// question b: bisection method
// [[Rcpp::export]]
double bisection(arma::dvec x, int a, int b){
  double theta;
  int numit = 0;
  while((b-a) > 1e-3 & numit < 10000){
  theta = (a+b)/2;
  numit = numit + 1;
  if (gprime1(x, a)*gprime1(x, theta) < 0) {
    b = theta;
  } else {
      a = theta;
      theta = (a+b)/2;}
  }
  return theta;
}


// question c: fixed-point iteration
// [[Rcpp::export]]
double fixedpoint(arma::dvec x, double init, double alpha){
  double theta = init;
  int numit = 0;
  int maxit = 100;
  while((std::fabs(gprime1(x, theta)) > 1e-8) & (numit < maxit)) {
  theta = theta + alpha*gprime1(x, theta);
  numit = numit + 1;
  }
  return theta;
}


// question d: secant method
// [[Rcpp::export]]
double secant(arma::dvec x, double theta0, double theta1){
  double theta = theta1;
  double theta_prev = theta0;
  int numit = 0;
  int maxit = 10000;
  while(std::fabs(gprime1(x, theta))<1e-8 | numit==maxit){
    theta = theta - gprime1(x, theta)*(theta-theta_prev)/(gprime1(x, theta)-gprime1(x, theta_prev));
    numit = numit + 1;
  }
  return theta;
}


// question e: Normal(theta, 1) (n=20)
// [[Rcpp::export]]
double n(arma::dvec x, const double theta){
  double value;
  value = -10*log(2*M_PI)-arma::sum(arma::pow((x-theta),2))/2;
  return value;
}

double nprime1(arma::dvec x, double theta){
  double value;
  value = arma::sum(x-theta);
  return value;
}

double nprime2(arma::dvec x, double theta){
  double value;
  value = -20.0;
  return value;
}

// [[Rcpp::export]]
double nnewton(arma::dvec x, double init, double eps){
  int i = 0;
  double theta = init;
  while(std::fabs(nprime1(x, theta)) > eps){
    theta = theta - nprime1(x, theta)/nprime2(x, theta);
    i = i + 1;
    // Rprintf("%.2d %.3f\n",i,theta);
  }
  return theta;
}

// [[Rcpp::export]]
double nbisection(arma::dvec x, int a, int b){
  double theta;
  int numit = 0;
  while((b-a) > 1e-3 & numit < 10000){
    theta = (a+b)/2;
    numit = numit + 1;
    if (nprime1(x, a)*nprime1(x, theta) < 0) {
      b = theta;
    } else {
      a = theta;
      theta = (a+b)/2;}
  }
  return theta;
}

// [[Rcpp::export]]
double nfixedpoint(arma::dvec x, double init, double alpha){
  double theta = init;
  int numit = 0;
  int maxit = 100;
  while((std::fabs(nprime1(x, theta)) > 1e-8) & (numit < maxit)) {
    theta = theta + alpha*nprime1(x, theta);
    numit = numit + 1;
  }
  return theta;
}

// [[Rcpp::export]]
double nsecant(arma::dvec x, double theta0, double theta1){
  double theta = theta1;
  double theta_prev = theta0;
  int numit = 0;
  int maxit = 10000;
  while(std::fabs(nprime1(x, theta))<1e-8 | numit==maxit){
    theta = theta - nprime1(x, theta)*(theta-theta_prev)/(nprime1(x, theta)-nprime1(x, theta_prev));
    numit = numit + 1;
  }
  return theta;
}


// question 2
// [[Rcpp::export]]
double f(arma::dvec x, double theta){
  double value;
  value = arma::sum(log(1-cos(x-theta)));
  return value;
}

double fprime1(arma::dvec x, double theta){
  double value;
  value = -arma::sum(sin(x-theta)/(1-cos(x-theta)));
  return value;
}

double fprime2(arma::dvec x, double theta){
  double value;
  value = -arma::sum(1/(1-cos(x-theta)));
  return value;
}

// [[Rcpp::export]]
double fnewton(arma::dvec x, double init, double eps){
  int i = 0;
  double theta = init;
  while(std::fabs(fprime1(x, theta)) > eps){
    theta = theta - fprime1(x, theta)/fprime2(x, theta);
    i = i + 1;
    // Rprintf("%.2d %.3f\n",i,theta);
  }
  return theta;
}



// additional hw: find maximum of f(x)=log(x)/(x+1) (1<x<5)
// [[Rcpp::export]]
double h(const double x){
  double value;
  value = std::log(x)/(x+1);
  return value;
}

double hprime1(const double x){
  double value;
  value = (1+(1/x)-std::log(x))/(std::pow((1+x),2));
  return value;
}

double hprime2(const double x){
  double value;
  value = (-1/(std::pow(x,2)+std::pow(x,3)))-2*(1+(1/x)-std::log(x))/(std::pow((1+x),3));
  return value;
}

// bisection method
// [[Rcpp::export]]
double hbisection(int a, int b){
  double x;
  int numit = 0;
  while((b-a) > 1e-3 & numit < 10000){
    x = (a+b)/2;
    numit = numit + 1;
    if (hprime1(a)*hprime1(x) < 0) {
      b = x;
    } else {
      a = x;
      x = (a+b)/2;}
  }
  return x;
}

// newton method
// [[Rcpp::export]]
double hnewton(double init, double eps){
  double x = init;
  int i = 0;
  double diff = 1.0;
  while(std::fabs(diff) > eps){
    diff= -hprime1(x)/hprime2(x);
    x = x + diff;
    i = i + 1;
    // Rprintf("%.2d %.3f\n",i,x);
  }
  return x;
}

// secant method
// [[Rcpp::export]]
double hsecant(double theta0, double theta1){
  double theta = theta1;
  double theta_prev = theta0;
  int numit = 0;
  int maxit = 10000;
  while(std::fabs(hprime1(theta))<1e-8 | numit==maxit){
    theta = theta - hprime1(theta)*(theta-theta_prev)/(hprime1(theta)-hprime1(theta_prev));
    numit = numit + 1;
  }
  return theta;
}

// fixed-point iteration
// [[Rcpp::export]]
double hfixedpoint(double init, double alpha){
  double theta = init;
  int numit = 0;
  int maxit = 100;
  while((std::fabs(hprime1(theta)) > 1e-8) & (numit < maxit)) {
    theta = theta + alpha*hprime1(theta);
    numit = numit + 1;
  }
  return theta;
}

